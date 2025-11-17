# INI file contents

```
[metadata]
Table Type = [Observations|]
Observation Type = [Spectroscopy|Speckle]
Instrument = [PEPSI|ZorroAlopeke|DSSI]
Analysis Type = [Speckle Detection|Data Assessment]

[transformations]
CoordinateDataTransformation = {ra_column}, {dec_column}
InferredTypeTransformation = {column_name}

[options]
constraint policy = [log][skip]

[columns]
{column_name} = str

[constraints]
{column_name} = {foreign_table}.{foreign_column}

[units]
{column_name} = {Astropy unit string}
```
