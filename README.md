# Writing analysis from scratch
Start here for setting up environment: [Example_of_using_DST_nodes (Wiki)](https://wiki.bnl.gov/sPHENIX/index.php/Example_of_using_DST_nodes#Building%20a%20package)
- Setup local compilation for bash shel:

```
source /opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=/sphenix/user/shulga/tpc2019_install
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL
```
- Creating package files:
`CreateSubsysRecoModule.pl --all --overwrite CalculateDistortions` 
(*very useful script providing all files needed for new package*)

- Compilation of an empty package:
```
mkdir build
cd build
../autogen.sh --prefix=$MYINSTALL
make -j 4
make install
```

Currently have a problem:
```
./.libs/libCalculateDistortions.so: undefined reference to `Fun4AllHistoManager::Fun4AllHistoManager(std::string const&)'
./.libs/libCalculateDistortions.so: undefined reference to `Fun4AllHistoManager::dumpHistos(std::string const&, std::string const&)'
```
- reading first file:
