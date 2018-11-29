

import rpy2.robjects as ro


import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()





codigo_r = """
heatmap_plot <- function(dataset) {
  library(pheatmap)
  pdf(file="Plot.pdf")
  profiling <- read.csv(dataset)
  m <- as.matrix(profiling[, -1])
  rownames(m) <- profiling$Comparison
  n.ind <- dim(profiling)[1]
  n_cluster<-dim(profiling$Comparison)[1]
  pheatmap(m, scale = "none", clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean", clustering_method = "average",
           cutree_rows = n.ind, fontsize = 11)
  dev.off()
}
"""
ro.r(codigo_r)
saluda_py = ro.globalenv['heatmap_plot']
res = saluda_py('file_to_plot.csv')
print(res[0])


codigo_r = """
saluda <- function(cadena) {
 return(paste("Hola, ", cadena))
}
"""
ro.r(codigo_r)

sentence = set()
if sentence.add("sff"):
    print("true")
if sentence.add("sffq"):
    print("true")
if sentence.add("sff"):
    print("true")
if sentence.add("sffq"):
    print("true")
if sentence.add("sffjkbd"):
    print("true")
if "sff" in sentence:
    print("esta")
for s in sentence:
    print(s)
value = ["?d1.", "qwpp", "sff", "sffq"]
sentence_where = ""
for l in value:
    if l not in sentence:
        sentence_where += l
    sentence.add(l)


if value[0].rstrip("."):
    print("djkdjkd")


if value[0].replace("1.", "2."):
    print(value[0])


ss = "qqqq89k0"
pp = "nmnm"

print(ss.lstrip("qqqq"))

print(ss.__contains__(pp))

dict = {}
dict["aa"] = "nmnm"
dict["as"] = "nmnm"
dict["ad"] = "nmnm"

if pp not in dict:
    print("OK")

#diccionario = {'color': 'rosa', 'marca': 'Zara', 'talle': 'U'}

diccionario = {"http://project-iasis.eu/vocab/gender": {
                "http://project-iasis.eu/Gender/C0025266": 0.6756756756756757,
                "http://project-iasis.eu/Gender/C0043210": 0.32432432432432434
            },
            "http://project-iasis.eu/vocab/hasDiagnosisTumorStage": {
                "http://project-iasis.eu/LCPatient/TumorStage/C0242379": 0.4594594594594595,
                "http://project-iasis.eu/LCPatient/TumorStage/C0001418": 0.24324324324324326,
                "http://project-iasis.eu/LCPatient/TumorStage/C0024121": 0.08108108108108109,
                "http://project-iasis.eu/LCPatient/TumorStage/C0152013": 0.13513513513513514,
                "http://project-iasis.eu/LCPatient/TumorStage/C0027651": 0.08108108108108109
            },
            "http://project-iasis.eu/vocab/tki": {
                "False": 1.0
            },
            "http://project-iasis.eu/vocab/hasDiagnosisAge": {
                "Age 66-75": 0.8378378378378378,
                "Age 56-65": 0.05405405405405406,
                "Age >75": 0.10810810810810811
            },
            "http://project-iasis.eu/vocab/brainMetastasis": {
                "False": 0.918918918918919,
                "True": 0.08108108108108109
            },
            "http://project-iasis.eu/vocab/familiarAntecedents": {
                "False": 0.6756756756756757,
                "True": 0.32432432432432434
            }}
claves = diccionario.keys()
print(claves)

array = [43, 435, 4359]
print(len(array))
print(array[len(array)-1])
if "http://project-iasis.eu/Gender/C0025266" in diccionario["http://project-iasis.eu/vocab/gender"]:
    print("esta")