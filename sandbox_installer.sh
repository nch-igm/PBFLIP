# Install mamba
conda install -n base -c conda-forge mamba
 
 # Change the directory to ~apps

 cd ~/apps

# install SQANTI3 and cDNA_Cupcake
pip install numpy pandas Cython

git clone https://github.com/Magdoll/cDNA_Cupcake.git
cd cDNA_Cupcake

export PATH=$PATH:/data/apps/cDNA_Cupcake/sequence/
export PATH=$PATH:/data/apps/cDNA_Cupcake/rarefaction/
python setup.py build
python setup.py install

chmod a+x /data/apps/cDNA_Cupcake/sequence/*.py

# Install SQANTI3
cd ..

git clone https://github.com/ConesaLab/SQANTI3.git

# Some SQANTI3 dependancies

mamba install -c bioconda BCBioGFF

R -e "install.packages('reshape', repos='http://cran.rstudio.com/')"
R -e "install.packages('gridExtra', repos='http://cran.rstudio.com/')"
R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')"
R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/');BiocManager::install('NOISeq')"
R -e "install.packages('ggplotify', repos='http://cran.rstudio.com/')"
R -e "install.packages('markdown', repos='http://cran.rstudio.com/')"
R -e "install.packages('DT', repos='http://cran.rstudio.com/')"
R -e "install.packages('plotly', repos='http://cran.rstudio.com/')"

# gtfToGenePred will be downloaded to ~/bin

wget -O ~/bin/gtfToGenePred https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
chmod a+x ~/bin/gtfToGenePred


export PATH=$PATH:$HOME/bin:$HOME/apps/SQANTI3:$HOMW/apps/cDNA_Cupcake/sequence/
export PYTHONPATH=$HOME/miniconda3/lib/python3.7/site-packages/:$HOME/apps/cDNA_Cupcake/sequence/:$HOME/apps/cDNA_Cupcake/cupcake/tofu/

wget https://github.com/broadinstitute/picard/releases/download/2.26.1/picard.jar -O $HOME/bin/picard.jar 

chmod a+x $HOME/bin/picard.jar

# SNPeff

# Go to home dir
cd $HOME/apps/

# Download latest version
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

# Unzip file
unzip snpEff_latest_core.zip
rm  snpEff_latest_core.zip 
cd snpEff

sed -ie 's/\.\/data\/\/data\/g' snpEff.config
java -jar snpEff.jar download GRCh38.99