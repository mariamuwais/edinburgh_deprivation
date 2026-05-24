# edinburgh_deprivation
# Exploring Deprivation in Rural Edinburgh
An analysis of deprivation patterns in Edinburgh using the Scottish Index of Multiple Deprivation indicators, looking at how the drivers of deprivation vary between urban and rural areas and proposing a modified deprivation index with differentiated domain weightings for urban and rural contexts.
## Method
### Geographically Weighted Regression (GWR)
GWR was applied to explore how relationships between individual SIMD domain indicators and overall deprivation rank vary spatially across Edinburgh. 
### Random Forest Classification and Regression
Random forest models were trained separately for urban and rural data zones to identify which SIMD domain indicators are most predictive of deprivation rank in each context. Variable importance scores were used to compare the relative contribution of each domain across the two area types.
### Modified Deprivation Index
Building on the findings from GWR and random forest analysis, a modified SIMD is proposed that applies differentiated domain weights for urban and rural areas, reflecting the differing ways in which deprivation is experienced and driven across settlement types.
## Data
Scottish Index of Multiple Deprivation (SIMD).
Data zone boundaries (https://www.data.gov.uk/dataset/afe1aca6-bea2-4283-8847-eecf98ab41a4/data-zone-boundaries-2022). 
Urban/Rural Classification (https://www.gov.scot/publications/scottish-government-urban-rural-classification-2022/).
## Author
Developed as part of an MSc course at the University of Edinburgh.
