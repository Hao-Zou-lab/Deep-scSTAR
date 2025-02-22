![Image text](https://github.com/Hao-Zou-lab/Deep-scSTAR/blob/main/logo.png)
# **Deep scSTAR documentation**

Deep Learning of single-cell State Transition Across-samples of RNA-seq data (**Deep scSTAR**) ,  is a deep learning model designed for mapping high-dimensional single-cell RNA-seq data into a lower-dimensional feature space, enabling the distinction of cell phenotypes and extraction of dynamic, phenotype-associated features while preserving the integrity of original biological signals.
![Image text](https://github.com/Hao-Zou-lab/Deep-scSTAR/blob/main/concept.png)
## **Dependency**

```shell
    python >= 3.10.0
    pytorch >= 2.2.1
    scikit-learn >= 1.2.2
    numpy >= 1.25.2
    pandas >= 1.5.3
    seaborn >- 0.13.1
    matplotlib >= 3.5.2
    imbalanced-learn >= 0.10.1
    umap-learn >= 0.5.5
```

We recommend three primary platforms for installing and running Deep scSTAR:

1. **Google Colab**: Given that Deep scSTAR was developed on Colab, running it on this platform allows you to bypass the environment setup process entirely. This method is highly recommended for users looking for quick deployment and execution without the need to manage dependencies and environments locally.
2. **Linux**: For users preferring local installation on Linux, it's advised to install dependencies in a sub-environment using miniconda3. 
3. **Windows via Anaconda Prompt**: For Windows users, utilizing Anaconda Prompt offers a straightforward method to run Deep scSTAR. After setting up the conda environment with all necessary dependencies, you can execute the Deep scSTAR scripts directly in the prompt.



## **Installation and Running Deep scSTAR**

### Running on Google Colab

Please check [Demo-Colab.ipynb](https://github.com/Hao-Zou-lab/Deep-scSTAR/blob/main/Demo_Colab.ipynb) for the usage of Deep scSTAR. The demo takes around 11 minutes to run on a server (GPU: T4 GPU). It demonstrates how to read `adata` using `scanpy`, generate inputs for Deep scSTAR, and run the analysis.

### Reproducing the Results from the Paper

For users interested in reproducing the results presented in our paper, we have provided two additional notebooks in the [experiments](https://github.com/Hao-Zou-lab/Deep-scSTAR/tree/main/experiments) directory:

1. **LPC-Responsive Endothelial Cell Data**  
   Notebook: [Code_to_Reproduce_Results_for_LPC_Responsive_Endothelial_Cell_Subpopulations.ipynb](https://github.com/Hao-Zou-lab/Deep-scSTAR/blob/main/experiments/Code_to_Reproduce_Results_for_LPC_Responsive_Endothelial_Cell_Subpopulations.ipynb)  
   **Description:** This code is used to reproduce the results for the LPC-Responsive Endothelial Cell dataset.  
   **Requirements:**  
   - The notebook needs to be run on Google Colab with an A100 GPU for optimal performance.  
   **Additional Resources:**  
   - At the end of the notebook, you will find links to pretrained models and the main results.

2. **B-cell Simulation Data**  
   Notebook: [Code_to_Reproduce_Results_for_the_B_cell_Simulation_Data.ipynb](https://github.com/Hao-Zou-lab/Deep-scSTAR/blob/main/experiments/Code_to_Reproduce_Results_for_the_B_cell_Simulation_Data.ipynb)  
   **Description:** This code is used to reproduce the results for the B-cell simulation dataset.  
   **Requirements:**  
   - The notebook needs to be run on Google Colab with an A100 GPU for optimal performance.  
   **Additional Resources:**  
   - At the end of the notebook, you will find links to pretrained models and the main results.

Simply request a Google Colab server with the required A100 GPU, run these notebooks, and you will obtain results consistent with those reported in the paper.

### Running on Linux via Miniconda3 Terminal

If you do not have Miniconda3, install it following the official guidelines. After installing Miniconda3, you can create a new environment, for instance, named DscSTAR (you can choose any name you prefer). Open your terminal and create an environment called DscSTAR:

```shell
  conda create -n DscSTAR python=3.10
```

Activate your environment, then retrieve the Deep scSTAR package from GitHub. If network issues prevent the use of `git clone`, you can download the  [DscSTAR.zip](https://github.com/Hao-Zou-lab/Deep-scSTAR/archive/refs/heads/main.zip) file directly and extract it to your local path. You can directly click on the [case](https://drive.google.com/file/d/1-XWWTZzxaw5GsIPUk4wPUshYg29QC-bk/view?usp=drive_link) and [ctr](https://drive.google.com/file/d/1N5qTWe-LIMmIiiobjKbtmm4Ju4ug0WYf/view?usp=drive_link) links to download the demo data, or use `gdown` to download the data:

```shell
  codeconda activate DscSTAR
  git clone https://github.com/Hao-Zou-lab/Deep-scSTAR.git
  # If you need to change the current drive, use the command "cd /d YOUR-PATH".
  cd Deep-scSTAR
  # download demo data
  pip install gdown
  
  gdown https://drive.google.com/uc?id=1pLvi62-PcIPMHAVlOlMym6xhWPZtRiKu -O Deep-scSTAR/inputs/case.h5ad
  gdown https://drive.google.com/uc?id=1-41wQpuQDpugk4sF0nQmfoQnaQYjOvql -O Deep-scSTAR/inputs/ctr.h5ad
```

Before setting up the environment, remember to check your CUDA version compatibility with the PyTorch version you plan to install for GPU acceleration. The `nvidia-smi` command will help you determine your CUDA version.

Set Up Environment:

```shell
  pip install matplotlib==3.7.1
  pip install scikit-learn==1.2.2
  # Check compatibility of CUDA version with PyTorch to use GPU acceleration
  pip install torch torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cu121
  pip install imbalanced-learn==0.10.1
  pip install pandas==1.5.3
  pip install seaborn==0.13.1
  pip install umap-learn==0.5.5
```

To run Deep scSTAR, execute the following command in your terminal:

```shell
python run_DscSTAR.py --f "h5ad" --s "h5ad" \ # --f sets input format, --s sets output format, both as H5AD.
--data-folder "inputs/" \
--input-1 "case.h5ad" \
--input-2 "ctr.h5ad" \
--output-1 "case.out.h5ad" \
--output-2 "ctr.out.h5ad" \
--output-folder "outputs/"
```

The demo takes around 5 minutes to run (GPU: A100).



### Running on Windows via Anaconda Prompt

If you do not have  Anaconda, install it. After installing Anaconda, you can create a new environment, for example, DscSTAR (you can change to any name you like). Open Anaconda Prompt create an environment called DscSTAR

```shell
  conda create -n DscSTAR python=3.10
```

Activate your environment, get the python package from github. If network issues prevent the use of `git clone`, you can directly download the  [DscSTAR.zip](https://github.com/Hao-Zou-lab/Deep-scSTAR/archive/refs/heads/main.zip) file and extract it locally. You can directly click on the [case](https://drive.google.com/file/d/1-XWWTZzxaw5GsIPUk4wPUshYg29QC-bk/view?usp=drive_link) and [ctr](https://drive.google.com/file/d/1N5qTWe-LIMmIiiobjKbtmm4Ju4ug0WYf/view?usp=drive_link) links to download the demo data, or use `gdown` to download the data:

```shell
  conda activate DscSTAR
  git clone https://github.com/Hao-Zou-lab/Deep-scSTAR.git
  # If you need to change the current drive, use the command "cd /d YOUR-PATH".
  cd ./DscSTAR
  # download demo data
  pip install gdown
  gdown https://drive.google.com/uc?id=1pLvi62-PcIPMHAVlOlMym6xhWPZtRiKu -O Deep-scSTAR/inputs/case.h5ad
  gdown https://drive.google.com/uc?id=1-41wQpuQDpugk4sF0nQmfoQnaQYjOvql -O Deep-scSTAR/inputs/ctr.h5ad
```

Before setting up the environment, it's important to note that when installing PyTorch, you must ensure that the version is compatible with your CUDA version to utilize GPU acceleration. If you're uncertain about your CUDA version, you can check it using the `nvidia-smi` command in the command prompt.

Set Up Environment:

```shell
  pip install matplotlib==3.7.1
  pip install scikit-learn==1.2.2
  pip3 install torch torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cu121
  pip install imbalanced-learn==0.10.1
  pip install pandas==1.5.3
  pip install seaborn==0.13.1
  pip install umap-learn==0.5.5
```

Finally, to run Deep scSTAR, use the following command:

```shell
python run_DscSTAR.py --f "h5ad" --s "h5ad" \ # --f sets input format, --s sets output format, both as H5AD.
--data-folder "inputs/" \
--input-1 "case.h5ad" \
--input-2 "ctr.h5ad" \
--output-1 "case.out.h5ad" \
--output-2 "ctr.out.h5ad" \
--output-folder "outputs/"
```

The demo takes around 6 minutes to run (GPU: RTX 4090).



### Outputs

After running `run_DscSTAR.py`, you will obtain five output files in the `outputs/` directory, which include:

1. `case.out.h5ad` and `ctr.out.h5ad`: These files contain the reconstructed data.
2. `UMAP.png`: This image presents a UMAP plot of the reconstructed data.
3. `Loss.png`: This graph displays the average total loss versus the number of training epoch.
4. `trained_model.pth`: This is the trained model file, which stores the learned weights and biases of the neural network after training.




## **Using .csv file as inputs**

If you wish to use .csv file as inputs, we provide a method to do so after reading the data of interest in R. You can add the phenotypic labels you desire. The following is a demonstration of generating demo data in R:

```R
HSP_expression <- FetchData(SeuratObj, vars = "HSP90AA1")
HSP_info <- ifelse(HSP_expression > 0, "HSP+", "HSP-")
SeuratObj[["HSPinfo"]] <- HSP_info

Ctr_filtered <- subset(seurat_filtered, subset = HSPinfo == 'HSP-')
Case_filtered <- subset(seurat_filtered, subset = HSPinfo == 'HSP+')

Ctr_filtered <- as.matrix(Ctr_filtered@assays$RNA@data)
Case_filtered <- as.matrix(Case_filtered@assays$RNA@data)

n <- ncol(Ctr_filtered)
colnames(Ctr_filtered) <- 1:n
colnames(Case_filtered) <- (n+1):(n + ncol(Case_filtered))

write.csv(Ctr_filtered, file="ctr.csv", row.names=TRUE)
write.csv(Case_filtered, file="case.csv", row.names=TRUE)
```

In this example, `SeuratObj` is a Seurat object containing normalized gene expression data. Let's define cells expressing the gene `HSP90AA1` as the case group, and those not expressing it as the control group. Then format the column names in a sequential manner and output the data to CSV files. The `case.csv` and `ctr.csv` files are structured to meet the input requirements of Deep scSTAR, with gene names as row headers and columns labeled with sequential numbers from `1` to `m` (where `m` is the total number of cells across both files). Columns in `ctr.csv` are numbered from `1` to `n`, and in `case.csv`, from `n+1` to `m` (with `n` being the number of columns in `ctr.csv` and `m` the total column count in both files, noting that `ctr` and `case` labels could be interchanged but the column naming from `1` to `m` is essential), ensuring each cell has a unique identifier across both files. This arrangement includes a top header row showing these numeric column identifiers and gene names in the leftmost column, aligning with Deep scSTAR's specifications for gene expression data analysis.

```
python run_DscSTAR.py --f "h5ad" --s "h5ad" \ # --f sets input format, --s sets output format, both as csv.
--data-folder "inputs/" \
--input-1 "case.csv" \
--input-2 "ctr.csv" \
--output-1 "case.out.csv" \
--output-2 "ctr.out.csv" \
--output-folder "outputs/"
```

The output `case.csv` and `ctr.csv` can be used as inputs for Deep scSTAR. Afterwards, you can read the processed data and merge it with the metadata of the original Seurat object for downstream analysis.

```R
file_path1 <- "./case.out.csv"
case <- as.matrix(read.csv(file_path1, header=TRUE, sep=",", row.names=1))
file_path2 <- "./ctr.out.csv"
ctr <- as.matrix(read.csv(file_path2, header=TRUE, sep=",", row.names=1))

colnames(case) <- colnames(Case_filtered)
colnames(ctr) <- colnames(Ctr_filtered)
combined_matrix <- cbind(case, ctr)

SeuratObj_processed <- CreateSeuratObject(counts = combined_matrix, project = "SeuratObj_processed")
SeuratObj_processed@meta.data <- SeuratObj@meta.data[match(rownames(SeuratObj_processed@meta.data), rownames(SeuratObj@meta.data)), ]
SeuratObj_processed[["RNA"]]@data <- SeuratObj_processed[["RNA"]]@counts
```




## **Adjustable Model Parameters**

The following parameters can be adjusted within `run_DscSTAR.py`:

`batch_size` - The number of samples per mini-batch. Select a number that is approximately 5% of the total number of samples.

`num_epochs` - The number of iterations for training the model. You can choose between 200-400. A higher number may lead to overfitting of the model.

`gamma` - If the results appear distorted, you can appropriately increase the gamma value. A higher gamma value can constrain the model to retain more of the original information.



## Applying pre-trained model

You can apply your model that has already been trained on a  new dataset by running apply_pretrained.py:

```shell
  python apply_pretrained.py --input-file1 "case.csv" --input-file2 "ctr.csv" \
            --output-name1 "case.applyed.csv" --output-name2 "ctr.applyed.csv" \
            --data-folder "./inputs/" \
            --model-path "./outputs/trained_model.pth"
```
After running `apply_pretrained.py`, you will obtain three output files in the `outputs/applyed_model/` directory, which include:

1. `case.applyed.csv` and `ctr.applyed.csv`: These files contain the reconstructed data.
2. `UMAP.png`: This image presents a UMAP plot of the reconstructed data.

