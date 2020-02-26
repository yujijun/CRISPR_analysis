args <- '~/Jobs/Data_Archive/CRISPR_Screens/ImmuneScreens/ImmuneScreens_Annotation.csv'
# args <- "/Users/yujijun/Dropbox (Partners HealthCare)/ImmuneScreens/ImmuneScreens_Annotation.csv"
anno_file <- args[1]
dirname <- dirname(anno_file)
annodata <- read.csv(anno_file, stringsAsFactors = FALSE)
annodata <- annodata[annodata$CountFile!="", ]

#### Construct the prefix for each study ####
prefix = paste0(annodata$Dir_name, ".", annodata$Model, "_", annodata$Condition)
idx = prefix %in% prefix[duplicated(prefix)]
tmp = paste0(prefix[idx], "_r") #mark all replicated file
for(i in unique(tmp)){ tmp[tmp==i] = paste0(i, 1:sum(tmp==i)) }
prefix[idx] = tmp; prefix = gsub("^\\w*\\.", "", prefix)
prefix <- gsub(" ","",prefix)
dirname <- file.path(dirname, annodata$Dir_name)
#### Design matrix ####
get_design <- function(anno){# c("Model", "Condition", "Day0", "Control", 'Treatment')
  anno = as.vector(anno)[c("Model", "Condition", "Day0", "Control", 'Treatment')]
  Day0 = unlist(strsplit(as.character(anno[3]), ", *"))
  Control = unlist(strsplit(as.character(anno[4]), ", *"))
  Treatment = unlist(strsplit(as.character(anno[5]), ", *"))
  if(all(Day0=="None")){
    tmp = data.frame(row.names = c(Control, Treatment),
                     cond = rep(0:1,c(length(Control), length(Treatment))))
    design = model.matrix(~cond, tmp)
    colnames(design) = c("baseline", paste0(anno[1], "_", anno[2]))
  }else{
    design = data.frame(row.names = c(Day0, Control, Treatment), Baseline = 1,
                     Control = rep(c(0,1,0), c(length(Day0), length(Control), length(Treatment))),
                     Treatment = rep(c(0,0,1), c(length(Day0), length(Control), length(Treatment))))
    colnames(design) = c("baseline", gsub(",.*", "", as.character(anno[4])),
                         gsub(",.*", "", as.character(anno[5])))
    colnames(design)[-1] = paste0(anno[1], "_", anno[2], colnames(design)[-1])
  }
  design = cbind(Sample = rownames(design), design)
  return(design)
}

# output all designmatrix
design_path = file.path(dirname, paste0("designmatrix/", prefix, "_design.txt"))
for (i in seq(dim(annodata)[1])){
  dir.create(dirname(design_path[i]), showWarnings = FALSE)
  design_matrix <- get_design(annodata[i,])
  write.table(design_matrix, file=design_path[i], sep='\t',
              row.names = FALSE, col.names = TRUE, quote = FALSE)
}

#### Run MAGeCK ####
mkdir_cmd = paste0("mkdir ", file.path(rep(dirname, each=3), c("count", "test", "mle")))
count_cmd = paste0("mageck count -k ", dirname, "/rawdata/", annodata$CountFile,
                   " -l ", dirname, "/lib/", annodata$library_calculation,
                   " --norm-method median -n ", dirname, "/count/",
                   annodata$Dir_name, "_", prefix)
test_cmd = paste0("mageck test -k ", dirname, "/rawdata/", annodata$CountFile,
                  " -c ", annodata$Control, " -t ", annodata$Treatment,
                  " --remove-zero-threshold 5 --sort-criteria pos",
                  " --norm-method median -n ", dirname, "/test/", prefix)
mle_cmd = paste0("mageck mle -k ", dirname, "/rawdata/", annodata$CountFile,
                 " -d ", design_path, " --remove-outliers --threads 8 ",
                 "--norm-method median -n ", dirname, "/mle/", prefix)

## Change the normalization method
idx = grepl("norm", annodata$CountFile)
count_cmd[idx] = gsub("norm-method median", "norm-method none", count_cmd[idx])
test_cmd[idx] = gsub("norm-method median", "norm-method none", test_cmd[idx])
mle_cmd[idx] = gsub("norm-method median", "norm-method none", mle_cmd[idx])
count_cmd = count_cmd[!duplicated(paste(annodata$Dir_name, annodata$CountFile))]
# system(paste(count_cmd,collapse = ";"))
# system(paste(test_cmd,collapse = ";"))
# system(paste(mle_cmd,collapse = ";"))

bash = c(count_cmd, test_cmd, mle_cmd)
write.table(bash, file.path(dirname(anno_file), "quick_runmageck.sh"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

