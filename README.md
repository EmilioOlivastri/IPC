# README #

This is the official repository for the paper *IPC: Incremental Probabilistic Consensus-based Consistent Set Maximization for SLAM Backends*.

## Prerequisites

- G2O (Tag [20201223_git](https://github.com/RainerKuemmerle/g2o/releases/tag/20201223_git))
- yaml-cpp
- Eigen3
- Boost
- Python3

## How to build

Remember to change the *G2O_ROOT* variable in the CMakeLists file to the correct position on your machine.

```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc --all)
```

## How to use

1. Select your favorite dataset pose graph optimization from [here](https://lucacarlone.mit.edu/datasets/). The dataset has to be in g2o format, some of them may need conversion.

2. Generate the ground truth file from that dataset using [SE-Sync](https://github.com/david-m-rosen/SE-Sync). Generate in such a way that the output is a text file that lists the poses of the trajectory in the following format (x, y, theta)

3. Spoil the dataset of your choice with the desired number of outliers using:
```
python3 scripts/generateDataset.py -i /path/to/clean/g2o_file -n n_outliers 
```
> > The original version of this script is from the [vertigo](https://github.com/OpenSLAM-org/openslam_vertigo/blob/master/datasets/generateDataset.py) package.

4. Adjust the config file based on the examples of the cfg folder.

```
canonic_inliers : 1614 <- Number of correct loop closures;
fast_reject_th : 6.251 <- Chi2 threshold for first check;
fast_reject_iter_base : 50 <- Number of optimization steps for first check;
slow_reject_th : 6.251 <- Chi2 threshold for second check;
slow_reject_iter_base : 100 <- Number of optimization steps for second check;
```
5. Run the following command for running the tester for IPC:
```
./ipc_tester -c ../cfg/INTEL_params.yaml
```
It will produce 2 files: output_name.txt and output_name.PR.
The first file contains the final estimated trajectory (x, y, theta), while the latter
contains (precision, recall, average time of convergence).

ROS was used for visualization of the obtained trajectories, while for the estimation of the ATE/RPE a modified version from [Hipe](https://github.com/rvp-group/srrg2-hipe/tree/main).

## Contact information

- Emilio Olivastri [emilio.olivastri@phd.unipd.it](mailto:emilio.olivastri@phd.unipd.it)

## Cite Us

if you find this implementation and/or research helpful, please consider to cite:
```bibtex
@INPROCEEDINGS{olivastri2024ipc,
  title={{IPC}: Incremental Probabilistic Consensus-based Consistent Set Maximization for SLAM Backends},
  author={Olivastri, Emilio and Pretto, Alberto},
  booktitle={2024 IEEE International Conference on Robotics and Automation (ICRA)},
  year={2024}
}
```

## License
Distributed under the BSD 2 License. See ```LICENSE``` for more information.
