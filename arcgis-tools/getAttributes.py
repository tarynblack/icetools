# ArcGIS script to fill glacier outline attributes
# T Black, 29 July 2020

import arcpy
import attributes as atb

null_value = str(-999)

# Set parameters. Outlines are chosen by user in ArcGIS GUI.
outlines = arcpy.GetParameterAsText(0)
python_version = 'PYTHON_9.3'
field_names = [f.name for f in arcpy.ListFields(outlines)]

# Create a cursor object to loop through and update each row in the table.
with arcpy.da.UpdateCursor(outlines, field_names) as cursor:
    # For each row, get the Image_ID and use it to update each other field.
    for row in cursor:
        if row[field_names.index('Image_ID')]:
            image_id = row[field_names.index('Image_ID')]

        if 'Sensor' in field_names and atb.getSensor(image_id) != null_value:
            row[field_names.index('Sensor')] = atb.getSensor(image_id)
        if 'Image_Tile_Coordinate' in field_names and \
            atb.getTileCoordinate(image_id) != null_value:
            row[field_names.index('Image_Tile_Coordinate')] = \
                atb.getTileCoordinate(image_id)
        if 'Image_Reference_System' in field_names and \
            atb.getReferenceSystem(image_id) != null_value:
            row[field_names.index('Image_Reference_System')] = \
                atb.getReferenceSystem(image_id)
        if 'Source_Date' in field_names and \
            atb.getSourceDate(image_id) != null_value:
            row[field_names.index('Source_Date')] = atb.getSourceDate(image_id)
        if 'Circadian_Date' in field_names and \
            atb.getCircadianDate(image_id) != null_value:
            row[field_names.index('Circadian_Date')] = \
                atb.getCircadianDate(image_id)
        if 'Year' in field_names and atb.getYear(image_id) != null_value:
            row[field_names.index('Year')] = atb.getYear(image_id)
        if 'Season' in field_names and atb.getSeason(image_id) != null_value:
            row[field_names.index('Season')] = atb.getSeason(image_id)
              
        cursor.updateRow(row)

del cursor, row
