mssm<-read.csv("C:\\Data\\lung.csv", header=T)
head(mssm)
mssm$LR_T<-ifelse(mssm$LR_T<0, mssm$TERM_T, mssm$LR_T)
# LR_T < 0: TERM_T, LR_T >= 0: LR_T
mssm$META_T<-ifelse(mssm$meta==1,mssm$META_T, mssm$TERM_T)
# meta(전이) == 1: MERA_T, meta == 0: TERM_T
mssm$death<-2-mssm$death
# death==1: 죽음, death==0:생존
mssm$LRMeta.s<-ifelse(mssm$meta+mssm$LR==2,1,0)
# 전이O & 국소재발O: 1, 아니면 0
mssm$LR_T<-ifelse(mssm$LRMeta.s==1&mssm$LR_T==mssm$META_T,mssm$LR_T-0.5, mssm$LR_T)
# 전이O & 국소재발O & LR_T=META_T: LR_T - 0.5
mssm$LRMeta_T<-ifelse(mssm$LRMeta.s==1, apply(cbind(mssm$LR_T,mssm$META_T),1,max),mssm$TERM_T)
# 전이O & 국소재발O: 각 행의 LR_T와 META_T 중 더 큰 값 
mssm$LR_T<-ifelse(mssm$LR_T==0&mssm$LR==1,0.5, mssm$LR_T)
# LR_T ==0 & 국소재발O: 0.5
mssm$META_T<-ifelse(mssm$META_T==0&mssm$meta==1,0.5, mssm$META_T)
# META_T==0 & 전이O: 0.5

install.packages("mstate")
library(mstate)

# 전이확률 행렬
tmat<-transMat(x=list(c(2,3,5), c(4),c(4,5),c(5), c()), names=c("Op","LR","Meta","LRMeta","Death"))
tmat

mssm1<-msprep(data=mssm, trans=tmat, time=c(NA,"LR_T","META_T","LRMeta_T","TERM_T"), status=c(NA,"LR","meta","LRMeta.s","death"),
              keep = c("sex","smoke","age","adjuv"))
mssm1[mssm1$Tstart==mssm1$Tstop,]$Tstop<-mssm1[mssm1$Tstart==mssm1$Tstop,]$Tstop+0.5 
mssm1[, c("Tstart","Tstop","time")]<-mssm1[, c("Tstart","Tstop","time")]/365.25
mssm1[mssm1$id < 3, ]
covs<-c("sex","smoke","age","adjuv")
mssm1<-expand.covs(mssm1,covs,longnames=FALSE,append=TRUE)


#########################################################################################################################################
##### model1: full model
c1 <- coxph(Surv(Tstart, Tstop, status) ~ sex.1 + sex.2 + sex.3 + sex.4 + sex.5 + sex.6 + sex.7
               + smoke.1 + smoke.2 + smoke.3 + smoke.4 + smoke.5 + smoke.6 + smoke.7
               + age.1 + age.2 + age.3 + age.4 + age.5 + age.6 + age.7
               + adjuv.1 + adjuv.2 + adjuv.3 + adjuv.4 + adjuv.5 + adjuv.6 + adjuv.7 + strata(trans), 
               data = mssm1, method = "breslow")
c1
c1$loglik[2]
1- pchisq(c1$loglik[2]-c1$loglik[2],2)


##### model2: 공변량을 축소한 모형
# full model에서 유의하지 않았던 공변량들을 하나의 변수로 합쳐보자
mssm2<- mssm1
mssm2$sex12367 <- ifelse(mssm2$trans %in% c(1,2,3,6,7), mssm2$sex, 0)
# mssm2$smoke는 만들 필요 X
mssm2$age12567 <- ifelse(mssm2$trans %in% c(1,2,5,6,7), mssm2$age, 0)
mssm2$adjuv13467 <- ifelse(mssm2$trans %in% c(1,3,4,6,7), mssm2$adjuv, 0)


# sex 공변량 축소
c11 <- coxph(Surv(Tstart, Tstop, status) ~ sex12367 +  sex.4 + sex.5
            + smoke.1 + smoke.2 + smoke.3 + smoke.4 + smoke.5 + smoke.6 + smoke.7
            + age.1 + age.2 + age.3 + age.4 + age.5 + age.6 + age.7
            + adjuv.1 + adjuv.2 + adjuv.3 + adjuv.4 + adjuv.5 + adjuv.6 + adjuv.7 + strata(trans), 
            data = mssm2, method = "breslow")
c11
1- pchisq(c1$loglik[2]-c11$loglik[2],2) # p-value가 크므로 간단한 모델이 더 적합함

# smoke 공변량 축소
c12 <- coxph(Surv(Tstart, Tstop, status) ~ sex.1 + sex.2 + sex.3 + sex.4 + sex.5 + sex.6 + sex.7
             + smoke
             + age.1 + age.2 + age.3 + age.4 + age.5 + age.6 + age.7
             + adjuv.1 + adjuv.2 + adjuv.3 + adjuv.4 + adjuv.5 + adjuv.6 + adjuv.7 + strata(trans), 
             data = mssm2, method = "breslow")

c12
1- pchisq(c1$loglik[2]-c12$loglik[2],2) # p-value가 작으므로 복잡한 모델이 더 적합함

# age 공변량 축소
c13 <- coxph(Surv(Tstart, Tstop, status) ~ sex.1 + sex.2 + sex.3 + sex.4 + sex.5 + sex.6 + sex.7
             + smoke.1 + smoke.2 + smoke.3 + smoke.4 + smoke.5 + smoke.6 + smoke.7
             + age12567 + age.3 + age.4 
             + adjuv.1 + adjuv.2 + adjuv.3 + adjuv.4 + adjuv.5 + adjuv.6 + adjuv.7 + strata(trans), 
             data = mssm2, method = "breslow")
c13
1- pchisq(c1$loglik[2]-c13$loglik[2],2) # p-value가 작으므로 복잡한 모델이 더 적합함

# adjuv 공변량 축소
c14 <- coxph(Surv(Tstart, Tstop, status) ~ sex.1 + sex.2 + sex.3 + sex.4 + sex.5 + sex.6 + sex.7
             + smoke.1 + smoke.2 + smoke.3 + smoke.4 + smoke.5 + smoke.6 + smoke.7
             + age.1 + age.2 + age.3 + age.4 + age.5 + age.6 + age.7
             + adjuv13467 + adjuv.2 + adjuv.5 + strata(trans), 
             data = mssm2, method = "breslow")
c14
1- pchisq(c1$loglik[2]-c14$loglik[2],2) # p-value가 크므로 간단한 모델이 더 적합함

# sex와 adjuv 공변량 둘 다 축소
c15 <- coxph(Surv(Tstart, Tstop, status) ~ sex12367 +  sex.4 + sex.5
             + smoke.1 + smoke.2 + smoke.3 + smoke.4 + smoke.5 + smoke.6 + smoke.7
             + age.1 + age.2 + age.3 + age.4 + age.5 + age.6 + age.7
             + adjuv13467 + adjuv.2 + adjuv.5 + strata(trans), 
             data = mssm2, method = "breslow")
c15
1- pchisq(c1$loglik[2]-c15$loglik[2],2) # p-value가 크므로 간단한 모델이 더 적합함


##### model3: 일부 전이가 같다고 가정한 모델
mssm3 <- mssm1
mssm3$strata1 <- mssm3$trans

# 전이 5 == 전이 4
mssm3$strata1[mssm3$trans==5] <- 4
# 모델 적합
c31 <- coxph(Surv(Tstart, Tstop, status) ~ sex.1 + sex.2 + sex.3 + sex.4 + sex.5 + sex.6 + sex.7
            + smoke.1 + smoke.2 + smoke.3 + smoke.4 + smoke.5 + smoke.6 + smoke.7
            + age.1 + age.2 + age.3 + age.4 + age.5 + age.6 + age.7
            + adjuv.1 + adjuv.2 + adjuv.3 + adjuv.4 + adjuv.5 + adjuv.6 + adjuv.7 + 
            + strata(strata1), data = mssm3, method = "breslow")
c31
c31$loglik[2]
c1$loglik[2]

# 전이 6 == 전이 7 == 전이 3
mssm3$strata2 <- mssm3$trans
mssm3$strata2[mssm3$trans==5] <- 4
mssm3$strata2[mssm3$trans %in% c(6,7)] <- 3
# 모델 적합
c32 <- coxph(Surv(Tstart, Tstop, status) ~ sex.1 + sex.2 + sex.3 + sex.4 + sex.5 + sex.6 + sex.7
             + smoke.1 + smoke.2 + smoke.3 + smoke.4 + smoke.5 + smoke.6 + smoke.7
             + age.1 + age.2 + age.3 + age.4 + age.5 + age.6 + age.7
             + adjuv.1 + adjuv.2 + adjuv.3 + adjuv.4 + adjuv.5 + adjuv.6 + adjuv.7 + 
               + strata(strata2), data = mssm3, method = "breslow")
c32
c32$loglik[2]
c1$loglik[2]

##### model4: intermediate events 고려한 모형
# 전이 6 == 전이 7 == 전이 3
mssm4 <- mssm1
mssm4$strata3 <- mssm3$trans
mssm4$strata3[mssm4$trans==5] <- 4
mssm4$strata3[mssm4$trans %in% c(6,7)] <- 3

# z2.6: 죽음 전 Meta 발생 여부
mssm4$z2 <- 0
mssm4$z2[mssm4$trans==6] <- 1
# z3.7: 죽음 전 LR+Meta 발생 여부
mssm4$z3 <- 0
mssm4$z3[mssm4$trans==7] <- 1
# z4.6: 죽음 전 Meta까지 걸린 시간
# z4.7: 죽음 전 LR+Meta까지 걸린 시간
mssm4$z4 <- 0
mssm4$z4[mssm4$trans %in% c(6,7)] <- mssm4$Tstart[mssm4$trans %in% c(6,7)]

covs<-c("z2","z3","z4")
mssm4<-expand.covs(mssm4,covs,longnames=FALSE,append=TRUE)

c4 <- coxph(Surv(Tstart, Tstop, status) ~ sex.1 + sex.2 + sex.3 + sex.4 + sex.5 + sex.6 + sex.7
             + smoke.1 + smoke.2 + smoke.3 + smoke.4 + smoke.5 + smoke.6 + smoke.7
             + age.1 + age.2 + age.3 + age.4 + age.5 + age.6 + age.7
             + adjuv.1 + adjuv.2 + adjuv.3 + adjuv.4 + adjuv.5 + adjuv.6 + adjuv.7 + 
             + z2.6 + z3.7 + z4.6 + z4.7 + strata(strata3), data = mssm4, method = "breslow")
c4
