import os
import re
import arcpy
import time
from functools import wraps



def time_now():
    t = time.localtime()
    current_time = time.strftime("%c", t)
    return current_time

def arcmsg27(func):
    '''
    Prints function arguments to console and arcgis tool messages. Use as decorator
    '''

    def time_now():
        t = time.localtime()
        current_time = time.strftime("%c", t)
        return current_time

    @wraps(func)
    def wrapper(*args, **kwargs):

        arcpy.AddMessage(''.join(['\n', time_now(), ' Start ', func.__name__]))
        for a in args:
            arcpy.AddMessage(''.join(['\t"', str(a), '"']))
        for kw in kwargs:
            arcpy.AddMessage(''.join(['\t', str(kw), ': "', str(kwargs[kw]), '"']))

        t1 = time.time()
        result = func(*args, **kwargs)
        seconds = time.time() - t1
        m, s = divmod(seconds, 60)
        h, m = divmod(m, 60)
        arcpy.AddMessage(
            ''.join(['Completed ', func.__name__, ' in: ', str(h), 'h ', str(m), 'm ', str(round(s, 2)), 's\n']))
        print('')

        return result

    return wrapper


# Add arcmsg27 functionality to arcpy tools
SpatialJoin_analysis = arcmsg27(arcpy.SpatialJoin_analysis)
Project_management = arcmsg27(arcpy.Project_management)
FeatureClassToFeatureClass_conversion = arcmsg27(arcpy.FeatureClassToFeatureClass_conversion)
Buffer_analysis = arcmsg27(arcpy.Buffer_analysis)

import arcpy
import os


def union_chain(layers, workspace, output_fc):
    garbage = []
    layer_index = list(enumerate(layers))
    union = layers[0]

    for i, layer in layer_index[1:-1]:  # skip last item. Final unioni will be output_fc

        print(layer)
        print('')
        new_union = os.path.join(workspace, f'Union_{os.path.basename(layer)}')
        arcpy.Union_analysis(
            [union, layer],
            new_union
        )

        union = new_union
        garbage.append(union)

        # Final union output
    arcpy.Union_analysis(
        [union, layers[-1]],
        output_fc
    )

    # Delete intermediate unions
    for fc in garbage:
        arcpy.Delete_management(fc)

    return output_fc

def has_featrues(feature_class):
    return int(arcpy.GetCount_management(feature_class).getOutput(0)) > 0

@arcmsg27
def select_copy(output_featureclass, input_featureclass, clip_features, layer_name, where, distance):
    """
    Selects input features intersecting clip_features and copies them to output FC
    """
    layer = arcpy.MakeFeatureLayer_management(
        in_features=input_featureclass,
        out_layer=layer_name,
        where_clause=where
    )

    arcpy.SelectLayerByLocation_management(
        in_layer=layer,
        overlap_type="WITHIN_A_DISTANCE",
        select_features=clip_features,
        selection_type="NEW_SELECTION",
        search_distance=distance
    )

    arcpy.CopyFeatures_management(
        in_features=layer,
        out_feature_class=output_featureclass
    )

    return output_featureclass

# Removed illegal chacters
def safe_string(text):
    """
    Makes arcgis friendly strings
    """
    text = re.sub(r',', r'_', text)
    return re.sub(r'[\s\\/-]|[\(\)]', r'', text)


class Constraint:
    """
    Holds data related to the constraint including intermediary data paths like projected and clipped data.
    """

    def __init__(self, name, setback, feature_class, where_clause, constraint_class):
        self.name = safe_string(name.title())
        self.setback = setback
        self.feature_class = feature_class
        self.where_clause = where_clause
        self.constraint_class = safe_string(constraint_class.title())
        self.clip_fc = None
        self.project_fc = None
        self.spatialjoin_fc = None
        self.constraint_fc = None
        self.constraint_shp = None
        self.constraint_setback_fc = None

    def attributes(self):
        """
        List of constraint attributes. Filters out intermediate data paths from object attributes
        """
        attr = {key: value for key, value in self.__dict__.items() if
                'shp' not in key and '_fc' not in key}
        return attr

    def _fields(self):
        return [key for key in self.__dict__.keys() if 'shp' not in key and '_fc' not in key]


class ConstraintsModel:
    """
    Processes constraints uses a list containing Constraint objects
    Creates intermediate gdbs for projected, clipped etc data
    """
    def __init__(self, sites_fc, constraints, output_folder, output_epsg=3400, overwrite_outputs=True):
        self.sites_fc = sites_fc
        self.constraints = constraints
        self.output_folder = output_folder
        self.output_epsg = output_epsg
        self.working_folder = os.path.join(self.output_folder, 'working')
        self.clip_gdb = os.path.join(self.working_folder, 'Clip.gdb')
        self.merge_gdb = os.path.join(self.working_folder, 'Merge.gdb')
        self.dissolve_gdb = os.path.join(self.working_folder, 'Dissolve.gdb')
        self.project_gdb = os.path.join(self.working_folder, 'Project.gdb')
        self.spatialjoin_gdb = os.path.join(self.working_folder, 'SpatialJoin.gdb')
        self.gdb_folder = os.path.join(self.output_folder, 'gdb')
        self.shp_folder = os.path.join(self.output_folder, 'shp')
        self.constraint_classes = set([c.constraint_class for c in self.constraints])
        self.overwrite_outputs = overwrite_outputs
        self._validate_constraint_names(self.clip_gdb)

    @arcmsg27
    def process_constraints(self):
        self._create_output_gdbs()
        for constraint in self.constraints:
            self.process_constraint(constraint)

    def process_constraint(self, constraint):
        """
        Selects and copies constraint features within setback distance of self.sites_fc
        Projects data to out_epsg
        Used spatial join to add SiteID to outputs
        Buffers constraint features = setback distance
        """
        arcpy.AddMessage(constraint.name)
        # Select and copy features within the setback distance of the project. Add 10m to the selection just in case
        # of transformation issues
        constraint.clip_fc = os.path.join(self.clip_gdb, constraint.name)

        self._check_overwrite_fc(
            select_copy,
            constraint.clip_fc,
            output_featureclass=constraint.clip_fc,
            input_featureclass=constraint.feature_class,
            clip_features=self.sites_fc,
            layer_name=constraint.name,
            where=constraint.where_clause,
            distance=f"{int(constraint.setback) + 10} Meters"
        )

        # If nothing selected/clipped skip further processing
        if has_featrues(constraint.clip_fc):
            constraint.project_fc = os.path.join(self.project_gdb, constraint.name)

            self._check_overwrite_fc(
                Project_management,
                constraint.project_fc,
                in_dataset=constraint.clip_fc,
                out_dataset=constraint.project_fc,
                out_coor_system=self.output_epsg
            )

            # Add Site ID (assumues no overlapping features)
            constraint.constraint_fc = os.path.join(
                self.gdb_folder,
                constraint.constraint_class + '.gdb',
                constraint.name
            )

            self._check_overwrite_fc(
                SpatialJoin_analysis,
                constraint.constraint_fc,
                target_features=constraint.project_fc,
                join_features=self.sites_fc,
                out_feature_class=constraint.constraint_fc,
                join_operation='JOIN_ONE_TO_ONE',
                join_type='KEEP_ALL',
                match_option='CLOSEST'
            )

            self._add_constraints_attributes(constraint)

            # buffer if there is a setback & export
            dissolve_fields = list(constraint.attributes().keys()) + ['SiteID']
            if constraint.setback > 0:
                constraint.constraint_setback_fc = os.path.join(
                    self.gdb_folder,
                    constraint.constraint_class + '.gdb',
                    constraint.name + f"_{str(constraint.setback).replace('.', '_')}m"
                )

                self._check_overwrite_fc(
                    arcpy.Buffer_analysis,
                    constraint.constraint_setback_fc,
                    in_features=constraint.constraint_fc,
                    out_feature_class=constraint.constraint_setback_fc,
                    buffer_distance_or_field=constraint.setback,
                    dissolve_option='LIST',
                    dissolve_field=dissolve_fields
                )

                constraint.dissolve_fc = os.path.join(
                    self.dissolve_gdb,
                    constraint.name + f"_{str(constraint.setback).replace('.', '_')}m"
                )
                self._check_overwrite_fc(
                    arcpy.Dissolve_management,
                    constraint.dissolve_fc,
                    in_features=constraint.constraint_setback_fc,
                    out_feature_class=constraint.dissolve_fc,
                    dissolve_field=dissolve_fields,
                    multi_part='MULTI_PART'
                )

            else:
                constraint.dissolve_fc = os.path.join(self.dissolve_gdb, constraint.name)
                self._check_overwrite_fc(
                    arcpy.Dissolve_management,
                    constraint.dissolve_fc,
                    in_features=constraint.constraint_fc,
                    out_feature_class=constraint.dissolve_fc,
                    dissolve_field=dissolve_fields,
                    multi_part='MULTI_PART'
                )
        else:
            print('NO FEATURES SKIPPING')

    @arcmsg27
    def export_shp(self):
        gdbs = [os.path.join(self.gdb_folder, f) for f in os.listdir(self.gdb_folder) if f.endswith('.gdb')]
        for gdb in gdbs:
            arcpy.env.workspace = gdb
            out_folder_name = os.path.basename(gdb).replace('.gdb', '')
            out_folder = os.path.join(self.shp_folder, out_folder_name)
            os.makedirs(out_folder, exist_ok=True)
            arcpy.FeatureClassToShapefile_conversion(arcpy.ListFeatureClasses(), out_folder)

    def _validate_constraint_names(self, workspace):
        invalid_names = []
        for constraint in self.constraints:
            new_name = arcpy.ValidateTableName(constraint.name, workspace)
            if constraint.name != new_name:
                invalid_names.append((constraint.name, new_name))
                constraint.name = new_name
        arcpy.AddWarning(f"Invalid constraint names updated {invalid_names}")

    def _create_output_gdbs(self):
        os.makedirs(self.working_folder, exist_ok=True)
        os.makedirs(self.gdb_folder, exist_ok=True)

        # Create working gdbs
        for gdb in [self.clip_gdb, self.project_gdb, self.spatialjoin_gdb, self.dissolve_gdb, self.merge_gdb]:
            if not os.path.exists(gdb):
                arcpy.CreateFileGDB_management(
                    os.path.dirname(gdb),
                    os.path.basename(gdb)
                )
        # Create output gdbs
        for c in self.constraint_classes:
            c_gdb = os.path.join(self.gdb_folder, c + '.gdb')
            if not os.path.exists(c_gdb):
                arcpy.CreateFileGDB_management(
                    os.path.dirname(c_gdb),
                    os.path.basename(c_gdb)
                )

    def _check_overwrite_fc(self, func, output_fc, *args, **kwargs):
        """
        Checks if func should execute based on self.overwrite_outputs
        """
        execute = True
        if not self.overwrite_outputs and arcpy.Exists(output_fc):
            execute = False
        if execute:
            result = func(*args, **kwargs)
        else:
            arcpy.AddMessage(''.join(["\n", time_now(), ' ', func.__name__, ' outputs already exist ', output_fc]))
            result = output_fc
        return result

    def _create_shp_folders(self):
        for c in self.constraint_classes:
            os.makedirs(os.path.join(self.shp_folder, c), exist_ok=True)

    def _add_constraints_attributes(self, constraint):
        for attr, value in constraint.attributes().items():
            arcpy.CalculateField_management(
                in_table=constraint.constraint_fc,
                field=attr,
                expression=f'r"{value}"',
                field_type='TEXT'
            )

    def union_constraints(self, ordered_constraints):

        arcpy.env.workspace = self.dissolve_gdb
        unioned_fc = os.path.join(self.merge_gdb, f"Union_ConstraintsAll")

        union_chain(
            arcpy.ListFeatureClasses(),
            self.dissolve_gdb,
            unioned_fc
        )

        unioned_fc = self._merge_union_attributes(unioned_fc, ordered_constraints)

        return unioned_fc

    def _merge_union_attributes(self, union_fc, ordered_constraints):

        invalid_constraint_classes = set(ordered_constraints) - self.constraint_classes
        if len(invalid_constraint_classes) > 0:
            raise ValueError(f"Invalid constraint class: {invalid_constraint_classes}")

        fields = []
        for constraint in self.constraints:
            fields += constraint._fields()
        merge_attributes = list(set(fields))

        union_fields = [f.name for f in arcpy.ListFields(union_fc)]
        for attribute in merge_attributes:
            concat_fields = [f for f in union_fields if attribute in f]

            expression_fields = "!" + "!, !".join(concat_fields) + "!"
            expression = f"', '.join(sorted([f for f in [{expression_fields}] if f]))"

            arcpy.CalculateField_management(
                in_table=union_fc,
                field=attribute,
                expression=expression,
                expression_type='PYTHON3',
                field_type="TEXT"
            )

        arcpy.Dissolve_management(
            in_features=union_fc,
            out_feature_class=union_fc + "_dissolved",
            dissolve_field=merge_attributes
        )

        c_classes = [f"'{c}'" for c in ordered_constraints]
        c_classes_string = ', '.join(c_classes)

        code = f"""def val(c):
            for c_class in [{c_classes_string}]:
                if c_class in c:
                    return c_class"""

        arcpy.CalculateField_management(
            in_table=union_fc + "_dissolved",
            field="ConstraintCategory",
            expression='val(!constraint_class!)',
            expression_type='PYTHON3',
            code_block=code,
            field_type="TEXT"
        )


        arcpy.Intersect_analysis([self.sites_fc, union_fc + "_dissolved"], union_fc + "_Sites")

        return union_fc + "_dissolved"
