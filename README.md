
# 🧬 Normalization of RNA-Seq data using adaptive trimmed mean with multi-reference

This project compares various normalization methods for RNA-seq differential expression (DE) analysis on **Sequencing quality control project (SEQC)**, **Microarray quality control project dataset (MAQC)**, **Pickrell dataset** (`SRP001540`) from the `recount` project and two Simulated datasets. Multiple normalization strategies are applied, and performance is evaluated using ROC curves, DE genes and False discovery rate (FDR).


Vikas Singh, Nikhil Kirtipal, Byeongsop Song, Sunjae Lee, Normalization of RNA-Seq data using adaptive trimmed mean with multi-reference, Briefings in Bioinformatics, Volume 25, Issue 3, May 2024, bbae241, https://doi.org/10.1093/bib/bbae241

## 🛠️ Installation

You can download the source code directly from GitHub:

🔽 **[Download ZIP](https://github.com/vikkyak/Normalization-of-RNA-seq-data/archive/refs/heads/main.zip)**

Or clone the repository using Git:

```bash
git clone https://github.com/vikkyak/Normalization-of-RNA-seq-data.git




📚 Requirements
recount, edgeR, DESeq2, ROC, pROC, PoissonSeq, TCC, compositions, psych
install.packages(c("pROC", "psych"))
BiocManager::install(c("recount", "edgeR", "DESeq2", "PoissonSeq", "TCC", "compositions"))

👨‍💻 Author
Vikas Singh
📧 vikkysingh07@gmail.com


🧬 Citation
If you use this code or method in your research, please cite:

**Vikas Singh, Nikhil Kirtipal, Byeongsop Song, Sunjae Lee**  
**Normalization of RNA-Seq data using adaptive trimmed mean with multi-reference**  
*Briefings in Bioinformatics*, Volume 25, Issue 3, May 2024, bbae241  
🔗 [https://doi.org/10.1093/bib/bbae241](https://doi.org/10.1093/bib/bbae241)


