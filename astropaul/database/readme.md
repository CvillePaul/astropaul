# INI file contents

```
[metadata]
Table Type = Observations
Observation Type = [Spectroscopy|Speckle]
Instrument = [PEPSI|ZorroAlopeke|DSSI]
Analysis Type = [Speckle Detection|Data Assessment]
Associated Table = {table name}

[transformations]
CoordinateDataTransformation = {ra column}, {dec column}
InferredTypeTransformation = {column name}

[options]
Primary Key = {key generator class}|[class options]
    SequentialIntKeyGenerator
    ColumnBasedKeyGenerator|Column 1,Column 2[:4],Column 3[6:8]
Constraint Policy = [log][skip]

[columns]
{column name} = str

[constraints]
{column name} = {foreign table}.{foreign column}

[units]
{column name} = {Astropy unit string}
```
