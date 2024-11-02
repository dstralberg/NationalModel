library(ggplot2)
varimp <- read.csv("H:/Shared drives/BAM_NationalModels4/NationalModels4.0/Manuscript/variable importance/BAM_VariableImportance.csv")
varclass <- read.csv("H:/Shared drives/BAM_NationalModels4/NationalModels4.0/Manuscript/variable importance/BAM_VariableListClasses.csv")
varimp <- merge(varimp,varclass,by.x="variable",by.y="var")

varimpsum <- aggregate(varimp$importance, by=list("BCR"=varimp$region, "varclass" =varimp$varclass), FUN = "sum")
varlist <- aggregate(varimpsum$x, by=list("varclass"=varimpsum$varclass), FUN="sum")
names(varlist)[2] <- "total"
varimpbcr <- aggregate(varimpsum$x, by=list("BCR"=varimpsum$BCR), FUN="sum")
names(varimpbcr)[2] <- "total"
varimpsum <- merge(varimpsum,varimpbcr)
varimpsum$prop <- varimpsum$x/varimpsum$total

write.csv(varlist, file="H:/Shared drives/BAM_NationalModels4/NationalModels4.0/Manuscript/variable importance/BAM_VariableList.csv", row.names=FALSE)
write.csv(varimpsum, file="H:/Shared drives/BAM_NationalModels4/NationalModels4.0/Manuscript/variable importance/BAM_VariableImportanceSummary.csv", row.names=FALSE)

#p <- ggplot(varimpsum, aes(varclass, x)) + geom_bar(stat="identity") + facet_wrap(~ BCR, ncol=3)
#p <- ggplot(varimpsum) + geom_bar( aes(x = varclass, y = x, fill = varclass), stat = "identity") + facet_wrap(~ BCR, ncol=3)
p <- ggplot() + theme_bw() + geom_bar(aes(x=BCR, y=prop, fill=varclass), data=varimpsum ,stat = "identity") + scale_fill_manual(values=c("cadetblue","coral","black","gray","beige","darkolivegreen4","darkolivegreen1"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.title =element_blank()) 
png(filename="H:/Shared drives/BAM_NationalModels4/NationalModels4.0/Manuscript/variable importance/BAM_VariableImportance.png")
print(p)
dev.off()

varimpspp<- aggregate(varimp$importance, by=list("BCR"=varimp$region, "var" =varimp$variable, "species"=varimp$id), FUN = "sum")
varimpbcr <- aggregate(varimpsum$x, by=list("BCR"=varimpsum$BCR), FUN="sum")
names(varimpbcr)[2] <- "total"
varimpspp2 <- merge(varimpspp,varimpbcr)
varimpspp2$prop <- varimpspp2$x/varimpspp2$total
varimpsppprop <- aggregate(varimpspp2$prop, by=list("var" =varimpspp2$var, "species"=varimpspp2$species), FUN=mean)
varimpprop <- aggregate(varimpspp2$prop, by=list("var" =varimpspp2$var), FUN=mean)
write.csv(varimpprop, file="H:/Shared drives/BAM_NationalModels4/NationalModels4.0/Manuscript/variable importance/BAM_VariableImportanceMeanProp.csv", row.names=FALSE)
write.csv(varimpsppprop, file="H:/Shared drives/BAM_NationalModels4/NationalModels4.0/Manuscript/variable importance/BAM_VariableImportanceMeanPropBySpecies.csv", row.names=FALSE)
