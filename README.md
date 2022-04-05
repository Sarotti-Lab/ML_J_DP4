# ML_JDP4 Package

## An integrated quantum mechanics-machine learning approach for ultra-fast NMR structural elucidation

Authors: Ariel M. Sarotti & María M. Zanardi 

Usage: `ml_jdp4` 

### Installation Requirements

**ML_JDP4.py** needs python 3.8 or later to work. You could install the module by console using:
`pip3 install ml-jdp4`

### User Guide

You need to create a folder containing the following files: 

      1. The gaussian outputs of the NMR and NBO calculations (all conformers for all isomers).
      
      2. The excel file containing the experimental data and the labels of each nucleus associated with each experimental value.
      
 ### Technical requirements
 
**1) The output files:** must be named following the next convention: number_*.log or .out, where number represent the isomer number and could be from 1 to N where N is the number of candidate structures under study. For instance:
 
       1_NewNatProd_c01.log (Conformer 1 for isomer 1 of the compound names NewNatProd)

       1_NewNatProd_c02.log (Conformer 2 for isomer 1 of the compound names NewNatProd)

       2_NewNatProd_c01.log (Conformer 1 for isomer 2 of the compound names NewNatProd)

*The script allows the use of outputs form Gaussian03, 09 and 16.*

**2) The input excel file:** The experimental data and the labels of the candidate structures must be placed in an excel file following the next rules. The excel file should be constituted by two sheets; one containing the data for the coupling constants (named ‘J’) and the other with the NMR chemical shifts (named ‘shifts’).

**“J” sheet:** the first column *“exp”* should contain the experimental <sup>3</sup>J coupling constants. In case of interchangeable values, they should be arranged “upside-down” (that is, the larger value first). The second column *“exchange”* serves to indicate **0** (not interchangeable value) or **1** (experimental data interchangeable with its following value). The third and fourth columns are intended to place the labels of the coupled protons. For cases when there are isomers with different labels, there should be two columns for each isomer as indicated below.  

![image](https://user-images.githubusercontent.com/101136961/161282945-682190b8-2f04-4e53-bcbd-7e54b5dd9908.png){width=50%}

**“shifts” sheet:** the first column *“nuclei”* should contain the identity of the atom ‘c or C’ for <sup>13</sup>C and ‘h or H’ for hydrogen atoms. The second column *“exp_data”* should contain the experimental chemical shifts. In case of interchangeable values, they should be arranged *“upside-down”* (that is, the larger value first). The third column *“exchange”* serves to indicate **0** (not interchangeable value) or **1** (experimental data interchangeable with its following value). The following columns are intended to place the labels of the nuclei associated to the corresponding chemical shift. If 2 or more values are added in that region, the isotropic shielding values will be averaged (as in the case of methyl groups or equivalent methylene groups). For cases when there are isomers with different labels, there should be three columns for each isomer as indicated in the Figure.

![image](https://user-images.githubusercontent.com/101136961/161283203-35f3f2df-e6a3-43d4-b8b4-87eb0c7bca18.png)


**3) The output excel file:** once the ML-J-DP4.py is executed, a filed named *‘Results_ML_J_DP4.xlsx’* is created in the same folder. The file contains five sheets:

**DP4 sheet:**  tthe DP4 probabilities are shown for each isomer considering the information of H, C and J individually or altogether. Although the high accuracy in the ML predictions, it must be emphasized that some environments might not be correctly reproduced leading to large unscaled errors that would affect the scaling procedure and the concomitant J-DP4 values. Hence, to avoid potential misassignments, a sheet containing the unscaled shifts and errors are printed. In this scenario, it is advisable to recompute J-DP4 after removing or revising the conflicting signal.

![figura_dp4](https://user-images.githubusercontent.com/101136961/161762018-4b82f429-bfba-4b7a-ae08-8a9d8041d3dc.JPG)

**Unscaled Chemical Shifts sheet:** this sheet will be labeled as *“Shifts_Unsc”*, it displays the experimental chemical shifts (<sup>13</sup>C and <sup>1</sup>H) and the coupling constants (<sup>3</sup>J), in that order, along with the predicted values for each isomer.

![figura_Unsc](https://user-images.githubusercontent.com/101136961/161763338-f418d39b-4bfb-4469-9429-b7a0d4166b7b.JPG)

**Scaled Chemical Shifts sheet:** this sheet will be labeled as *“Shifts”*, it displays the experimental chemical shifts (<sup>13</sup>C and <sup>1</sup>H) and the coupling constants (<sup>3</sup>J), in that order, along with the predicted values after the scaling process for each isomer.

![figura_shifts](https://user-images.githubusercontent.com/101136961/161764360-20ec494b-f4ea-4dae-a8c8-5d566fbd710c.JPG)

**Unscaled Errors sheet:** this sheet will be labeled as *“Unsc_Errors”*, it displays the differences (absolute value) between Unscaled and experimental chemical shifts and coupling constants. It is necessary to check this sheet to verify that there are no nucleus with excessively large errors that could lead to an incorrect assignment.

![figura_UnsError](https://user-images.githubusercontent.com/101136961/161769513-056d7800-9052-4aa0-9de6-cc021beafc3f.JPG)

**Errors sheet:** it displays the differences (absolute value) between scaled and experimental chemical shifts and coupling constants. 

![figura_errors](https://user-images.githubusercontent.com/101136961/161770080-f0033838-f1ac-4459-a64f-e19aecd08dc7.JPG)

# Case study: (-)-Menthol

In order to illustrate the ML-dJ-DP4 and the ML-iJ/dJ-DP4 workflows, we present the analysis of (-)-menthol following these two approaches. As indicated in the Figure 1, there are four possible isomers. 
 
 ![mentol_1](https://user-images.githubusercontent.com/101136961/161816781-85c9528c-3053-447b-b773-f7eed058f5d9.JPG)

**Figure 1**

Following the recommended computational procedure, a total number of 123 conformers were found at the MMFF force field. Each structure was submitted to NMR calcualtions at the RHF/STO-3G level (with the pop=nbo option). The corresponding output files are provided in the Folder **“menthol_ML_dJ-DP4”**. According to Gaussian numbering scheme, the labels corresponding to each nuclei is given in Figure 2. 

 ![mentol_2](https://user-images.githubusercontent.com/101136961/161816816-7a39a084-1014-4d9d-8465-478f3b0511f0.JPG)

**Figure 2.** Carbon label followed by its corresponding proton(s) label(s) between parenthesis

## ML-dJ-DP4

The input excel file required for running the ML-J-DP4.py script is filled out as follows. The Excel file is also provided in the folder *“menthol_ML_dJ-DP4”*.
**J sheet**
 
![image](https://user-images.githubusercontent.com/101136961/161816911-d9af40d2-7839-4dd3-b0e7-ef55f1a99e50.png)

**Shifts sheet**
 
 ![image](https://user-images.githubusercontent.com/101136961/161816946-ee6c2269-1daa-4103-a8b0-eb5f83940871.png)

Once the script is run, the resulting excel report file **“Results_ML_J_DP4”** will be generated. The output Excel file is provided in the folder **“menthol_ML_dJ-DP4”**.

**DP4 sheet:**

![image](https://user-images.githubusercontent.com/101136961/161816993-3d57ca86-2d45-4348-9387-96f2e3c5fd32.png)

## ML-iJ/dJ-DP4

This approach involves removing unsuitable conformations that are incompatible with selected experimental <sup>3</sup>J<sub>HH</sub> values. In this case, we selected the J<sub>16,19</sub>= 10 Hz to constrain the sampling, leading to 34 conformations.

![mentol_3](https://user-images.githubusercontent.com/101136961/161817451-1b4a3759-b336-480e-b6f7-e3ec8d4485a4.JPG)
 
**Figure 3.** J<sub>16,19</sub> (10 Hz) was chosen to constrain the conformations

Once removed the unsuitable conformations, the remaining conformations are submitted to GIAO NMR calculations at the RHF/STO-3G level with the pop=nbo option. The resulting output files are given in the folder **“menthol_ML_iJ-dJ-DP4”**. It is important to note that in this method, the selection of the conformations should be done before running the `ML-J-DP4.py script`. The Excel file with the experimental data is the same as in the dJ-DP4 analysis, but leaving out the information of the <sup>3</sup>J<sub>HH</sub> used to filter the conformations.

![image](https://user-images.githubusercontent.com/101136961/161817542-537ea267-c75d-4b58-81f1-498b7ccdf436.png)

**J sheet:**
 
Once the script is run, the resulting excel report file **“Results_ML_J_DP4” will be generated**. The output Excel file is provided in the folder *“menthol_ML_iJ-dJ-DP4”*.

![image](https://user-images.githubusercontent.com/101136961/161817581-b71c7b3e-d64f-491d-aa15-265886e58677.png)

**DP4 sheet:**
 


