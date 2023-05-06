exprMatr <- matrix(sample(1:10, 100, replace=TRUE), nrow=20)#在1~10中有放回的随机抽取100个数，分成20行创建矩阵
rownames(exprMatr) <- paste("Gene", 1:20, sep="")#给行命名
colnames(exprMatr) <- paste("Sample", 1:5, sep="")#给列命名
head(exprMatr)


library(GENIE3)
set.seed(123) #设置随机数种子，以便重复实验
weightMat <- GENIE3(exprMatr)
dim(weightMat)
weightMat[1:5,1:5]

regulators <- c(2, 4, 7)
regulators <- c("Gene2", "Gene4", "Gene7")
weightMat <- GENIE3(exprMatr, regulators=regulators) #指定基因使用GENIE3分析

regulatorsList <- list("Gene1"=rownames(exprMatr)[1:10],
                       "Gene2"=rownames(exprMatr)[10:20],
                       "Gene20"=rownames(exprMatr)[15:20]) #创建列表并命名
set.seed(123)
weightList <- GENIE3(exprMatr, nCores=1, targets=names(regulatorsList), regulators=regulatorsList, returnMatrix=FALSE)

weightMat <- GENIE3(exprMatr, treeMethod="ET", K=7, nTrees=50) #使用ET算法，选择7个regulators，生成50个树

set.seed(123)
weightMat <- GENIE3(exprMatr, nCores=4, verbose=TRUE)#指定4个核

linkList <- getLinkList(weightMat)
dim(linkList) #查看列表维度
head(linkList)

linkList <- getLinkList(weightMat, reportMax=5)#选择weight最大的5个

linkList <- getLinkList(weightMat, threshold=0.1)#筛选wiight>0.1，目标基因



