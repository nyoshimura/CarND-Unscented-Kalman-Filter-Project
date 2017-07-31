[//]: # (Image References)

[image1]: ./img/result.png "result"

# Unscented Kalman Filter Project
#### Self-Driving Car Engineer Nanodegree Program

In this project utilize an Unscented Kalman Filter to estimate the state of a moving object of interest with noisy lidar and radar measurements.
![alt text][image1]

Passing the project requires obtaining RMSE values that are lower that the tolerance outlined in the project [rubric](https://review.udacity.com/#!/rubrics/783/view).

There are 7 criterias:
---

1. My code should compile. (done)
2. px, py, vx, vy output coordinates must have an RMSE <= [.09, .10, 0.40, 0.30] (done)
3. My Sensor Fusion algorithm follows the general processing flow (done)
4. My Kalman Filter algorithm handles the first measurements appropriately. (done)
5. My Kalman Filter algorithm first predicts then updates. (done)
6. My Kalman Filter can handle radar and lidar measurements. (done)
7. My algorithm should avoid unnecessary calculations. (done)
