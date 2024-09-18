## About
Uses ESRI Arcpy\
Tools to read a list of constraint parameters from csv and perform processing on constraint features
* Select and copy input data near the area of interest
* Project inputs to common coordinate system
* Buffer features to create setbacks
* Union constraint layers
* Handles export to .shp and intermediate data locations
## Example
[Example Usage with Jupyter Notebook notebook](ExampleProcess.ipynb)
* Load constraint parameters from a csv including data path, setback distance, and constraint query
  * [Example input parameters in csv](Example_Constraints_Processing_Parameters.csv)
* Load parameters into Constraint and ConstraintModel classes for processing
* Execute processing and export to shp
