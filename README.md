# FloodModelDownscaling
An emulation approach inspired by Bayesian Reification but better than the first one



We're comparing downscaled, calibrated 30m projections to calibrated 10m projections. 
We also fit the model: flood_10m ~ poly(downscaled_flood_30m,2). 
We also see if any hydrological factors explain why some residuals are higher than others from this model.
