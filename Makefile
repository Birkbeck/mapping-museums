#
# This scraper is based on the scrapy tutorial:
#    http://doc.scrapy.org/en/latest/intro/tutorial.html
#
# Start URL: https://www.travelblog.org/North-America/United-States/California/Los-Angeles/
# install on Ubuntu with
#   $ sudo pip install scrapy
#

r: wait
	Rscript -e "rmarkdown::render('code/mapping_museums.Rmd',output_file = '../tmp/mapping_museums_r_out.html')" > tmp/r_log.txt 2>&1

conv:
	# sudo apt-get install gnumeric
	cd tmp/topicmodel && find . -maxdepth 1 -name "*.tsv" -exec ssconvert {} --export-type=Gnumeric_Excel:xlsx2 \;
	#ssconvert PyCrawl/tmp/ethnicelebscom_dataset.v7.gender.tsv  PyCrawl/tmp/ethnicelebscom_dataset.v7.gender.xlsx
	#ssconvert ethnicelebscom_dataset/ethnicelebscom_url_pro_ethno.tsv  ethnicelebscom_dataset/ethnicelebscom_url_pro_ethno.xlsx
	#ssconvert PyCrawl/tmp/ethnicelebscom_people_single.tsv  PyCrawl/tmp/ethnicelebscom_people_single.xlsx

apportion:
	Rscript code/mm_apportionment.R

backup:
	cd .. && zip -r -q -X ~andreaballatore/Dropbox/DRBX_Temp/052-MappingMuseums.backup.zip 052-MappingMuseums -x "*/.*"

clean_plots:
	-rm -r mm_ab_shared/plots;

share: backup clean_plots
	rsync -r --exclude=".*" plots mm_ab_shared

conv_img:
	TODO

clean:
	-@rm WebCrawlPy/*pyc WebCrawlPy/__pycache__/*pyc

wait:
	sleep 3;