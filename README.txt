
--------------------------------------------------------------------------------
Geant4 Compatibility
--------------------------------------------------------------------------------

Developed for Geant4.9.6


--------------------------------------------------------------------------------
Installation
-------------------------------------------------------------------------------- 

(replace /path/to/workdir with directory containing EicRichGemTB source directory)

mkdir /path/to/workdir/EicRichGemTB-build
mkdir /path/to/workdir/EicRichGemTB-install

cd /path/to/workdir/EicRichGemTB-build

cmake -DCMAKE_INSTALL_PREFIX=/path/to/workdir/EicRichGemTB-install /path/to/workdir/EicRichGemTB

make -jN (N = number of processor cores available)

make install -jN


--------------------------------------------------------------------------------
Running
-------------------------------------------------------------------------------- 

./richtb

