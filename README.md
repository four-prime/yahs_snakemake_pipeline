## Usage

The pipeline expects the following file structure. The `samples` directory should contain a single subdirectory for each sample to be run through the pipeline. Within each of these subdirectories, a single `.bam` file and a single `.fasta` or `.fa` file is expected. The base names of these files do not matter.  

**Note:** Gzipped files are not supported at this time.

├── config.json
├── README.md
├── samples
│   ├── some_sample
│   │   └── input
│   │       ├── some_file_a.bam
│   │       └── some_other_file_a.fasta
│   └── some_other_sample
│       └── input
│           ├── some_file_b.bam
│           └── some_other_file_b.fasta
├── downloads
└── Snakefile

---

Before running the pipeline, update the `sample_names` field in `config.json` to include the samples you wish to process.  
For the example structure above, this would look like:

```json
{
  "sample_names": [
    "some_sample",
    "some_other_sample"
  ]
}
```

---

To use Snakemake's SLURM job scheduling and containerization features, you need to create a Snakemake profile at:

`~/.config/snakemake/`

A template configuration file can be created using the following command (adjust paths and options as needed):

```bash
echo "executor: slurm
jobs: 100
cores: 32
latency-wait: 30
software-deployment-method: apptainer
apptainer-args: \"--bind /your/directory/of/intrest\"
default-resources:
  slurm_account: unix-slurm
  slurm_extra: \"'--auks=yes'\"
  partition: cpu
  runtime: 96h" > ~/.config/snakemake/slurm-apptainer/config.yaml
```

Refer to the Snakemake documentation for more details on profile configuration:
[Snakemake Configuration Documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html)


Before using snakemake's SLURM integration, install the required plugin:

```bash
pip install snakemake-executor-plugin-slurm
```


Finally, to start the pipeline, run the following command from the project root:


```bash
snakemake --profile slurm-apptainer
```

Once complete, your project file tree should look like so:

├── config.json
├── README.md
├── samples
│   ├── some_sample
│   │   ├── input
│   │   │   ├── some_file_a.bam
│   │   │   └── some_other_file_a.fasta
│   │   ├── work
│   │   └── output
│   └── some_other_sample
│       ├── input
│       │   ├── some_file_b.bam
│       │   └── some_other_file_b.fasta
│       ├── work
│       └── output
├── downloads
└── Snakefile

Within the output folder should be two files. A .fa file containing the final scaffold output of yahs, and a .hic file which can be passed to [Juicebox](https://github.com/aidenlab/Juicebox/wiki/Download) to visulaise Hi-C contact matrices for the assembly and perform manual corrections. Details on how to perform manual corrections will be provided at a future date.
