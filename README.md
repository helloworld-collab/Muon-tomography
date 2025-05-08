
## Contents

- **`lyon.cc`**  
  C++ script to process the ROOT file `DAT000003lyon.root` and generate histograms for muon `cos squared(θ)` and other parameters at various observation levels.

- **`chyulu.cc`**  
  Similar to `lyon.cc`, but processes the data from `DAT000003chyulu.root` corresponding to Chyulu Hills.

- **`DAT000003lyon.root`**  
  ROOT file containing simulated muon event data for Lyon. Produced by running CORSIKAPLOTTER and converting the output with a COAST-based parser.

- **`DAT000003chyulu.root`**  
  ROOT file containing muon event data for Chyulu Hills.

- **`corsika_plotter.cc`**  
  This is the core C++ program based on the COAST template (by Tanguy Pierog & Ralf Ulrich, 2020) modified to:
  - Read `.dat` files from CORSIKA
  - Extract muon data per observation level
  - Store information such as `Theta`, `Phi`, `ID`, and `Nmu` into ROOT TTrees and TNtuples
  - Save outputs as `.root` files for further plotting

## How to Use

1. Compile any of the `.cc` plotting files with ROOT:
   ```bash
   root -l lyon.cc
