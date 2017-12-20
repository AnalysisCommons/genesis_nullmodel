SHELL=/bin/bash
APP_NAME=genesis_nullmodel
ARCHIVE_DIR=$(APP_NAME)/resources/tmp
R_DIR=$(APP_NAME)/resources/home/dnanexus
SETUP_DIR=$(APP_NAME)/resources/home/dnanexus
MD5_ARCHIVE = setup/local/archives.md5

R_REPO_FILES = $(notdir $(wildcard R/*.R))
SETUP_REPO_FILES = $(notdir $(wildcard setup/cloud/*))

R_FILES = $(addprefix $(R_DIR)/, $(R_REPO_FILES))
SETUP_FILES = $(addprefix $(SETUP_DIR)/, $(SETUP_REPO_FILES))
DNA_NEXUS_FILES = $(addprefix $(APP_NAME)/, src/code.sh Readme.developer.md Readme.md dxapp.json )

.PHONY : all
all: $(SETUP_FILES) $(R_FILES) $(DNA_NEXUS_FILES) archives
 
$(APP_NAME): 
	mkdir -p $@/src
	mkdir -p $@/resources/tmp
	mkdir -p $@/resources/home/dnanexus


$(APP_NAME)/src/code.sh : app/code.sh $(APP_NAME)
	cp $< $@	


$(APP_NAME)/dxapp.json : app/dxapp.json $(APP_NAME)
	cp $< $@


$(APP_NAME)/Readme.developer.md : app/Readme.developer.md $(APP_NAME)
	cp $< $@


$(APP_NAME)/Readme.md : app/Readme.md $(APP_NAME)
	cp $< $@


$(R_DIR)/%.R :  R/%.R $(APP_NAME)
	cp $< $@


$(SETUP_DIR)/% : setup/cloud/% $(APP_NAME)
	cp $< $@
	
# Remove everything but archives 
.PHONY: clean 
clean:
	rm $(R_FILE) $(SETUP_FILES) $(DNA_NEXUS_FILES)
         

.PHONY : archives
archives: $(MD5_ARCHIVE) setup/local/check_update_archives.sh $(APP_NAME) 
	./setup/local/check_update_archives.sh $< $(ARCHIVE_DIR)


.PHONY : update_archive_md5sum 
update_archive_md5sum: $(ARCHIVE_DIR)
	md5sum $(ARCHIVE_DIR)/*.tar.gz | awk '{sub(".*/","AUX_DIRNAME/", $$2);print }' > $(MD5_ARCHIVE)


