# NetClass1.0.0: A noncontrol-dependent approach for microbiome subnetwork classification and central bacterial identification
## Introduction
`NetClass` is a novel microbial network classification framework to classify the key subnetwork and identify central microbes at any body site using a random walk algorithm and rank-sum ratio-entropy weight evaluation model. <br>
The workflow of NetClass is shown in the following Figure: (Details about this model can be found this paper.)
## News
- Nov, 2023: NetClass version 1.0.0 is launched.
## Installation
- Prerequisites: NetClass is developed under R (*version >= 3.6.1*).
- Latest version: The latest developmental version of NetClass can be downloaded from GitHub and installed from source by `devtools::install_github('YihuaWWW/NetClass1.0.0')`
## Manual
In the R terminal, please use the command `?NetClass` to access the help documents.
## Data Requirements
Two types of data are needed to run the NetClass:

- The relative abundance of microbial species in different individuals: tax_L7.txt
- The bipartite network made of correspondence reactions of species(or genus) towards phylum: phylum_species.txt 

## How to cite `NetClass`
Please cite the following manuscript:
> *NetClass: A noncontrol-dependent approach for microbiome subnetwork classification and central bacterial identification(2023).<br>* *Yihua Wang, Qingzhen Hou, Fulan Wei, Bingqiang Liu and Qiang Feng*
## License
- NetClass is licensed under the GNU General Public License v3.0.
- Improvements and new features of NetClass will be updated on a regular basis.
- If you have any questions or feedback, please contact Yihua Wang yihuawww@163.com.
