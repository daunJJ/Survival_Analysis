ventICU <- read.table("C:\\Data\\ventICU.txt", header=T)

####### 데이터탐색
head(ventICU)
length(unique(ventICU$id))
dim(ventICU)
table(ventICU$id)
table(ventICU$from, ventICU$to)  
unique(ventICU$enum)
ventICU[ventICU$id==1008422,]

## state 1, 2 사이 전이가 일어난 행만 출력
cycle_idx <- as.numeric(rownames(ventICU[(ventICU$from ==1 & ventICU$to == 2) | (ventICU$from == 2 & ventICU$to == 1), ]))
## 순환 전이 count
id_counts <- table(ventICU[cycle_idx,]$id)
max(id_counts)
sum(id_counts == 1)
sum(id_counts == 2)
sum(id_counts == 3)
sum(id_counts == 4)
sum(id_counts == 5)

## 중도절단 처리
ventICU$to<-ifelse(ventICU$to==999,"cens",ventICU$to)


######## 모델 적용
### 1. etm 패키지 사용
install.packages("etm")
library(etm)

ventICU2 <- ventICU

# time 변수는 estop으로 설정 
ventICU2$time <- ventICU2$estop
# tra(가능한 전이 행렬) 만들기
tra<-matrix(ncol=4, nrow=4,FALSE)
tra[1,2:4]<- TRUE
tra[2,c(1,3:4)] <-TRUE

ICU.etm <- etm(ventICU2,c("1","2","3","4"),tra, "cens", s= 0)

plot(ICU.etm)
library(lattice)
xyplot(ICU.etm)

clfs <- trprob(ICU.etm, "1 1") + trprob(ICU.etm, "1 4")
var.clfs <- trcov(ICU.etm, "1 1") + trcov(ICU.etm, "1 4") +
  + 2 * trcov(ICU.etm, c("1 1", "1 4"))

### 2. 순환 전이가 없는 8가지 상태로 변환
# 8가지 상태로 변환하기
ventICU3 <- ventICU
ventICU3$enum<-as.numeric(ventICU3$enum)
ventICU3$to<-as.numeric(ventICU3$to)
ventICU3$from<-as.numeric(ventICU3$from)

n=length(unique(ventICU3$id))
uni_id<-unique(ventICU3$id)
n_row<-nrow(ventICU3)

# 변환된 상태는 to2, from2로 표현
ventICU3$from2<-ventICU3$from
ventICU3$to2<-ifelse(ventICU3$to==3,7,ventICU3$to)
ventICU3$to2<-ifelse(ventICU3$to==4,8,ventICU3$to2)

# ntt는 고유 id별로 최대 enum을 구하기 위함
ntt=rep(0,n)

for (i in 1:n){
  for (j in 1:n_row){
    if(uni_id[i]==ventICU3$id[j]) ntt[i]=ventICU3$enum[j]
  }}

as.numeric(ntt)

# 순환 state를 개별 state로 변경
for (i in 1:(n_row-1)){
  if(is.na(ventICU3$to[i])==FALSE){  # censoring 되지 않은 경우
    ii=i+1
    if(ventICU3$enum[i]==1&ventICU3$from2[i]==2&ventICU3$to2[i]==1){
      # 전이 1번, ON에서 시작, OFF에서 끝
      ventICU3$to2[i]=3; ii=i+1; ventICU3$from2[ii]<-3
      }
    if(ventICU3$enum[i]>=2&ventICU3$to[i]<=2){  # 전이 2번 이상, ON 또는 OFF로 전이
      ventICU3$to2[i]<-ventICU3$from2[i]+1  # 순환이 아닌 다음 상태로 만들기
      ii=i+1
      if(ventICU3$id[i]==ventICU3$id[ii]){  # 전이가 한번 더 일어난 경우
        ventICU3$from2[ii]<-ventICU3$to2[i]  # 그때의 from state를 직전의 to state와 동일하게 변경
}}}}

table(ventICU3$to2)
table(ventICU3$from2)
table(ventICU3$from2, ventICU3$to2) 

ventICU3[ventICU3$id==1008422,]

# coxph 모형 적합
library(survival)
## 1.중도절단 처리가 가능하면 ...
# ventICU3$time <- ventICU3$estop - ventICU3$estart # time 변수 생성
sum(is.na(ventICU3$to2))  # 중도절단된 데이터 14개
#ventICU3$event <- ifelse(is.na(ventICU3$to2), 0, 1) # 중도절단되면 0 아니면 1
#ventICU3$to2<-factor(ventICU3$to2,levels=c(1,2,3,4,5,6,7,8))
#coxph_fit = coxph(Surv(time,event,to2)~age+sex.female, id=id, data = ventICU3)

## 2. 중도절단 처리가 불가능하면...
## 2-1. 중도절단은 999로 처리 -> 하나의 상태가 됨
ventICU3$to3 <- ifelse(is.na(ventICU3$to2), 999, ventICU3$to2)
ventICU3$to3<-factor(ventICU3$to3,levels=c(999,1,2,3,4,5,6,7,8)) # factor형으로 변경
sum(is.na(ventICU3$to3)) 
table(ventICU3$to3)
coxph_fit = coxph(Surv(estart, estop,to3)~age+sex.female, id=id, data = ventICU3)
summary(coxph_fit)
concordance(coxph_fit)$concordance # c-index: 
## 2-2. 중도절단 데이터 없애기
ventICU3$to2<-factor(ventICU3$to2,levels=c(1,2,3,4,5,6,7,8)) # factor형으로 변경
sum(is.na(ventICU3$to2)) 
table(ventICU3$to2)
coxph_fit2 = coxph(Surv(estart, estop,to2)~age+sex.female, id=id, data = ventICU3)
summary(coxph_fit2)
concordance(coxph_fit2)$concordance # c-index: 

# mstate 패키지를 사용하여 모델 적용
install.packages("mstate")
library(mstate)

# 전이확률행렬
tmat<-transMat(x=list(c(2,7,8), c(3,7,8),c(4,7,8),c(5,7,8), c(6,7), c(8), c(), c()), names=c("1","2","3","4","5","6","7","8"))
tmat

ventICU3$to2<-as.factor(ventICU3$to2)
mstate <- msprep(ventICU3, id = "id", time = "estop", status = "to2")

### 3. ventiliation on and off를 공변량으로 한 competing risk 분석
ventICU4 <- ventICU
ventICU4$cause<-ifelse(ventICU4$to==3,1,ifelse(ventICU4$to==4,2,0)) 
  # state1,2이면 0, state3이면 1, state4이면 2
table(ventICU4$cause)
head(ventICU4)

ventICU4$estart<-as.numeric(ventICU4$estart)
ventICU4$estop<-as.numeric(ventICU4$estop)

install.packages("cmprsk")
library(cmprsk)
cov<-ventICU4[,c(2:3,5)]  # 공변량 지정

tt<-as.numeric(ventICU4[,8])  # stop time
delta=ventICU4$cause  # 상태

crr(tt,delta,cov, failcode=1)
crr(tt,delta,cov, failcode=2)

# 그래프 그리기
install.packages("etm")
library(etm)
cif1 <- etmCIF(survival::Surv(estart, estop, cause!= 0) ~from, ventICU4,
               etype = cause, failcode =1)  # discharge
plot(cif1, ci.type = "bars", 
     col = c(1, 2), lty = 1, curvlab = c("Off", "On"))

cif2 <- etmCIF(survival::Surv(estart, estop, cause!= 0) ~from, ventICU4,
               etype = cause, failcode =2)  # death 
plot(cif2, ci.type = "bars", 
     col = c(1, 2), lty = 1, curvlab = c("Off", "On"))













