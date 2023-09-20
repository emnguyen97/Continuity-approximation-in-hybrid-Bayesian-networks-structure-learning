
#dat.dag = B.W2

#dat.dag = W2.f

## For partition MCMC
{
blklist= data.frame(from=rep(c("GW","BWZ", "BMI"),each=2),to=c("BM1","BM2"))
tem=data.frame(from=intersect(c("AC","SX","BMI","TV","EG","EM","CD","FS","FH","RP1","DP1","SL"),colnames(dat.dag)),to="BWZ")
blklist=rbind(blklist,tem)
tem=data.frame(from=intersect(c("SX","BMI","SE","AC", "TV","EG","EM","CD","INC", "FS","FH","RP1","DP1","SL","OD","GW","BWZ"),colnames(dat.dag)),to="FE1")
blklist=rbind(blklist,tem)

tem=data.frame(from=intersect(c("SX","BMI","SE","AC", "TV","EG","EM","CD","INC", "FS","FH","RP1","DP1","SL","OD","GW","BWZ"),colnames(dat.dag)),to="ME1")
blklist=rbind(blklist,tem)
tem=data.frame(from=intersect(c("SX","BMI","AC", "TV","EG","EM","CD", "FS","FH","RP1","DP1","SL","OD","BWZ"),colnames(dat.dag)),to="GW")
blklist=rbind(blklist,tem)

tem=data.frame(from=intersect(c("BMI","AC","SE", "TV","EG","EM","CD","INC","ME1","FE1","BM1","BM2", "FS","FH","RP1","DP1","SL","OD","GW","BWZ"),colnames(dat.dag)),to="SX")
blklist=rbind(blklist,tem)
tem=data.frame(from=intersect(c("BMI","AC","TV","EG","EM","CD","SL","OD","GW","BWZ"),colnames(dat.dag)),to="SE")
blklist=rbind(blklist,tem)
tem=data.frame(from=intersect(c("BMI"),colnames(dat.dag)),to="AC")
blklist=rbind(blklist,tem)

}


## Additional links for CLG model - structure MCMC
{
temp=data.frame(from=intersect(c("INC","SE","BM1","BM2","BWZ"),colnames(dat.dag)),to="AC")
blklist_clg=rbind(blklist,temp)
temp=data.frame(from=intersect(c("BM1","BM2"),colnames(dat.dag)),to="FE1")
blklist_clg=rbind(blklist_clg,temp)
temp=data.frame(from=intersect(c("BM1","BM2"),colnames(dat.dag)),to="ME1")
blklist_clg=rbind(blklist_clg,temp)
temp=data.frame(from=intersect(c("INC","SE","BM1","BM2"),colnames(dat.dag)),to="GW")
blklist_clg=rbind(blklist_clg,temp)
temp=data.frame(from=intersect(c("BMI","INC","SE","BM1","BM2","BWZ"),colnames(dat.dag)),to="DP1")
blklist_clg=rbind(blklist_clg,temp)
temp=data.frame(from=intersect(c("BMI","INC","SE","BM1","BM2","BWZ"),colnames(dat.dag)),to="FH")
blklist_clg=rbind(blklist_clg,temp)
temp=data.frame(from=intersect(c("BMI","INC","SE","BM1","BM2","BWZ"),colnames(dat.dag)),to="FS")
blklist_clg=rbind(blklist_clg,temp)
temp=data.frame(from=intersect(c("BMI","INC","SE","BM1","BM2","BWZ"),colnames(dat.dag)),to="ME2")
blklist_clg=rbind(blklist_clg,temp)
temp=data.frame(from=intersect(c("BMI","INC","SE","BM1","BM2","BWZ"),colnames(dat.dag)),to="OD")
blklist_clg=rbind(blklist_clg,temp)
temp=data.frame(from=intersect(c("BMI","INC","SE","BM1","BM2","BWZ"),colnames(dat.dag)),to="RP1")
blklist_clg=rbind(blklist_clg,temp)
temp=data.frame(from=intersect(c("BMI","INC","SE","BM1","BM2","BWZ"),colnames(dat.dag)),to="SL")
blklist_clg=rbind(blklist_clg,temp)
temp=data.frame(from=intersect(c("BMI","INC","SE","BM1","BM2","BWZ"),colnames(dat.dag)),to="FV")
blklist_clg=rbind(blklist_clg,temp)
temp=data.frame(from=intersect(c("BMI","INC","SE","BM1","BM2","BWZ"),colnames(dat.dag)),to="HF")
blklist_clg=rbind(blklist_clg,temp)
temp=data.frame(from=intersect(c("BMI","INC","SE","BM1","BM2","BWZ"),colnames(dat.dag)),to="HSD")
blklist_clg=rbind(blklist_clg,temp)
}

{
blklist= blackpartition(blklist,dat.dag)
blklist_clg= blackpartition(blklist_clg,dat.dag)}
