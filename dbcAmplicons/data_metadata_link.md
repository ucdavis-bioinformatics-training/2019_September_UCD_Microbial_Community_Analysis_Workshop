Experiment Metadata
===============================================

This document assumes [dbcAmplicons installing software](./dbcAmplicons_installing_software.md) has been completed.

**IF** for some reason it didn't finish, is corrupted or you missed the session, you can use my instance. Create or edit the file

	/share/workshop/$USER/mca_example/src/dbcA_profile

edit the contents (use nano) to use my environment.

<div class="script">export PATH=/share/workshop/msettles/mca_example/bin:$PATH
module load java/jdk1.8
export RDP_PATH=/share/workshop/msettles/mca_example/src/RDPTools
module load anaconda2
. /software/anaconda2/4.5.12/lssc0-linux/etc/profile.d/conda.sh
conda activate /share/workshop/msettles/mca_example/src/dbcA_virtualenv
export PYTHON_EGG_CACHE=/share/workshop/$USER/mca_example/src
</div>

Lets login and navigate to our directory

	cd /share/workshop/$USER/mca_example/

**1\.** Lets copy the workshop data into our workshop directory.

	cp -r /share/biocore/workshops/2019_Sept_MCA/Illumina_Reads /share/workshop/$USER/mca_example/.

Take a look at the files ... what is inside the Illumina_Reads folder?

Lets look at the first few reads of each file, below is an example line to view the first file.

	zless Illumina_Reads/MCAWorkshopDataset_S1_L001_R1_001.fastq.gz | head

Look at the head of the other files.

**2\.** Lets copy the workshop metadata into our workshop directory.

Next lets make a metadata directory and transfer our barcode, primer and sample sheet to their

	cd /share/workshop/$USER/mca_example
	cp -r /share/biocore/workshops/2019_Sept_MCA/metadata /share/workshop/$USER/mca_example/.


Take a look at the files ... what is inside the metadata folder?

View the metadata files in the folder

If all this is correct, we are ready to begin. If our environment isn't already sourced, lets do so and then test.

	source /share/workshop/$USER/mca_example/src/dbcA_profile
	dbcVersionReport.sh
