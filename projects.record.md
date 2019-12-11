## Project name： CRISPR SCREEN pipeline
## Time：05-12
## Project content：    
**1. update and sort out data and code(upload all data into kraken)**
**2. package all code for a input and a output**    

## Project execution process：
**0. upload all related code and data into kraken:**  
**1. Create superlink in this folder:**
```
/liulab/cwan/TiSig/static
```
**2. add variable into :**
```
(TiSig/src/immunesig/__init__.py)
```
```
BASE_DIR = os.path.dirname(os.path.realpath(__file__))
STATIC_DIR = os.path.join(BASE_DIR, '..', '..', 'static/‘) #
SCICB_DIR = os.path.join(STATIC_DIR, 'ScRNA', 'ICB’)   #
SCICB_RMD = os.path.join(BASE_DIR, 'scICB_evaluator', 'workflow.R’) #this is for R code 
```

**3. library(opt parse) to transfer parameter outside:**
```
Library(optparse)
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-s", "--species"), type="character", default="mm10",
              help="input species [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")\
,
  make_option(c("-d", "--data"), type="character", default="./",
              help="output file name [default= %default]", metavar="character")
);

#check input parameters
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

setwd(opt$data)
genelist <- opt$file
species <- opt$species
output <- paste0(opt$out, "/scICB_evaluator/")
scRNAPlot(genelist, species, output)
```
**4. add my module name into( src/immuneSig.py):**
```
parser.add_argument('--evaluators',nargs='+',default=['NonICB','ICB','Syngeneic','Screen', 'scICB'],help="Evaluators you want to run.")
```
**5. add load code into api:**
```
if 'scICB' in evaluators: 

logger.info( 

'Start evaluating association of input geneset with single-cell ICB data.') 

from immunesig import SCICB_RMD, SCICB_DIR

inputGeneList=os.path.join(request['output'],'input_geneset.txt') 

outputDire=request['output'] 

r_cmd = ' '.join(['Rscript',SCICB_RMD, '-f ',inputGeneList,'-s',species,'-o',request['output'], '-d', SCICB_DIR]) 

print(r_cmd) 

process=subprocess.Popen(r_cmd, shell=True).wait() 

print(process)#[0] is stdout

logger.info('Successfully finish single-cell ICB evaluator.') 
```
2.6 create a new branch and push
All done

## Results:  

## Conclusion:
