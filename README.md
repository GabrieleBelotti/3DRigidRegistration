# 3D Rigid Registration

Utility optimized for rigid registration between a CT and a CBCT.
CT is referred and is expected to be the Fixed Image. On the opposite side CBCT is expected to be the Moving Image

## Usage from command line

```bash
ImageRegistration3DClass <MetricNumber> CTName&Path CBCTName&Path <options>
```
Various options are provided and can be print out using the help command:

```bash
ImageRegistration3DClass -h
```

Suggested options are:
- Fixed cropping, activated with `-c` --> will crop the fixed image (CT) to the region that overlaps to the moving image (CBCT);
- Isocenter alignment, activated with `-rt <RTplanFilename>` --> will align the moving image to the first isocenter in the RTPlan file and execute option -c (this is under the assumption that the isocenter is a good initialization for this registration);
- `-res`, `-fres` and `-mres` --> call a resampling routine for one or both images, affecting accuracy and execution time.

## Notes

Currently the registration uses a BSpline interpolator for the registration process and also for the optional resamplings. Such choice slows down both processes but guarantee a higher degree of convergence.

At this stage the class is written to output and write and image. This could be optional, while the tranform should be fed to external routines

## Current developing

Multistage registration is currently implementable but kept locked as a non activated option. The goal would be to optimize the registration to give an acceptable improvement in short time. Creating a user varying initialization could be a good idea (select coarseness, number of steps and guarantee a continuity).

Feeding the output to a reinitialization routine.
