# notice: exe and output is in build file
# processing_cityjson
processing CityJSON calculate the volume, area, rough index and so on
processed file use cjio to convert them into correct format, so no need to worry lod/format issues
## empty data processing:
### skip it, if some building is empty we skip it
### return alternative data, if some data such as triangules after cdt is not exist, we return ring instead for surface orientation judgement(there is one case,where the ring exists but the triangules is not exist because two points overlab together so there is a silver)
### randomly assign, if orientation is on boundaries, we randomly assign it
### avoid zero, add extreme small value to avoid variable divided by zero
### omit it, some cases cgal help us handle it and just return null value/throw error(what ever, the program not crash), so we just skip it and return zero value for bad data
### assume it is right, since we dont need to correct incorrect data. we assume it is correct
### simplify it, no adjacent values for buildigns
## validation
except orientation, all index is the same as true value tested on cube three dimension is all 1
<img width="445" alt="image" src="https://user-images.githubusercontent.com/60132725/227620483-23303012-a288-4490-867e-3e5fc918923b.png">
orientation is check by check the direction using google map:
<img width="694" alt="image" src="https://user-images.githubusercontent.com/60132725/227620681-5f0aa69d-6cbb-448b-b963-c7fd4c137211.png">
<img width="868" alt="image" src="https://user-images.githubusercontent.com/60132725/227621227-dafdcd40-0c8e-45e6-9ade-2c31d5fa3259.png">
the best sample may be bk building
