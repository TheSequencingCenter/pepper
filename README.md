## P.E.P.P.E.R.
[![Build Status](https://travis-ci.com/kishwarshafin/pepper.svg?branch=master)](https://travis-ci.com/kishwarshafin/pepper)

`PEPPER` is a genome inference module based on recurrent neural networks that enables long-read variant calling and nanopore assembly polishing in the [PEPPER](https://github.com/kishwarshafin/pepper)-[Margin](https://github.com/UCSC-nanopore-cgl/margin)-[DeepVariant](https://github.com/google/deepvariant) pipeline. This pipeline enables nanopore-based variant calling with [DeepVariant](https://github.com/google/deepvariant).

<p align="center">
<img src="./img/PMDV_variant_calling_ONT_v5.png" alt="PEPPER-Margin-DeepVariant Variant Calling Workflow" width="720p"></img>
</p>

---
### Version 0.5 update

We have updated PEPPER to support `Guppy v5.0.7`. The current v0.5 only supports ONT variant calling.

**NOTE: The models of v0.5 are trained with Guppy v5.0.7 "sup" mode. If you use have Guppy v4.X data, [please use PEPPER v0.4](https://github.com/kishwarshafin/pepper/tree/r0.4)**

**If you want to use PEPPER-Margin-DeepVariant for assembly polishing or PacBio HiFi variant calling, [please use PEPPER v0.4](https://github.com/kishwarshafin/pepper/tree/r0.4)**

---

### How to cite
Please cite the following manuscript if you are using `PEPPER-Margin-DeepVariant`:


<details>
<summary><a href="https://www.biorxiv.org/content/10.1101/2021.03.04.433952v1"><b>bioRxiv:</b> Haplotype-aware variant calling enables high accuracy in nanopore long-reads using deep neural networks.</a></summary>
Authors: Kishwar Shafin, Trevor Pesout, Pi-Chuan Chang, Maria Nattestad, Alexey Kolesnikov, Sidharth Goel, <br/> Gunjan Baid, Jordan M Eizenga, Karen H Miga, Paolo Carnevali, Miten Jain, Andrew Carroll, Benedict Paten.
</details>

---
### Quickstart
Please follow the quickstart guides to assess your setup. Please follow case-study documentations for detailed instructions.
* **Docker**: [Oxford Nanopore variant calling quick start](./docs/quickstart/variant_calling_docker_quickstart.md).
* **Singularity**: [Oxford Nanopore variant calling quick start](./docs/quickstart/variant_calling_singularity_quickstart.md).

### Case studies

The variant calling pipeline can be run on [Docker](https://docs.docker.com/install/linux/docker-ce/ubuntu/) or [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html#quick-installation-steps). The case studies are designed on `chr20` of `HG002` sample.

The case-studies include input data and benchmarking of the run:
* Nanopore variant calling using **Docker**: [Link](./docs/pipeline_docker/ONT_variant_calling.md)
* Nanopore variant calling using **Singularity**: [Link](./docs/pipeline_singularity/ONT_variant_calling_singularity.md)
* Nanopore variant calling using **NVIDIA-docker**: [Link](./docs/pipeline_docker_gpu/ONT_variant_calling_gpu.md)

### License
[PEPPER license](./LICENSE), [Margin License](https://github.com/UCSC-nanopore-cgl/margin/blob/master/LICENSE.txt) and [DeepVariant License](https://github.com/google/deepvariant/blob/r1.1/LICENSE) extend to the trained models (PEPPER, Margin and DeepVariant) and container environment (Docker and Singularity).

### Acknowledgement
We are thankful to the developers of these packages:
* [htslib & samtools](http://www.htslib.org/)
* [pytorch](https://pytorch.org/)
* [ONNX](https://onnx.ai/)
* [hdf5 python (h5py)](https://www.h5py.org/)

### Authors
[PEPPER](https://github.com/kishwarshafin/pepper)-[Margin](https://github.com/UCSC-nanopore-cgl/margin)-[DeepVariant](https://github.com/google/deepvariant) pipeline is developed in a collaboration between [UC Santa Cruz genomics institute](https://ucscgenomics.soe.ucsc.edu/) and the [Genomics team in Google Health](https://health.google/health-research/genomics/).


### Fun Fact
<img src="https://vignette.wikia.nocookie.net/marveldatabase/images/7/72/Anthony_Stark_%28Earth-616%29_from_Iron_Man_Vol_5_2_002.jpg/revision/latest?cb=20130407031815" alt="guppy235" width="240p"> <br/>

The name "P.E.P.P.E.R." is inspired from an A.I. created by Tony Stark in the  Marvel Comics (Earth-616).

PEPPER is named after Tony Stark's then friend and the CEO of Resilient, Pepper Potts.
