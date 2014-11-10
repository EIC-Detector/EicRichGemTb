
--------------------------------------------------------------------------------
Geant4 Compatibility
--------------------------------------------------------------------------------

Developed for Geant4.9.6


--------------------------------------------------------------------------------
Installation
--------------------------------------------------------------------------------

Getting the code from GitHub

git clone https://github.com/nfeege/EicRichGemTb.git


--------------------------------------------------------------------------------
Installation
-------------------------------------------------------------------------------- 

(replace /path/to/workdir with directory containing EicRichGemTb source directory)

mkdir /path/to/workdir/EicRichGemTb-build
mkdir /path/to/workdir/EicRichGemTb-install

cd /path/to/workdir/EicRichGemTb-build

cmake -DCMAKE_INSTALL_PREFIX=/path/to/workdir/EicRichGemTb-install /path/to/workdir/EicRichGemTb

make -jN (N = number of processor cores available)

make install -jN (optional- you can run it from the build directory)


--------------------------------------------------------------------------------
Running
-------------------------------------------------------------------------------- 

./richtb

