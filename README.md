<h1 align="center"> 🧪 Survival Analysis </h1>
<h4 align="center"> 의생명 통계 방학 학부인턴 과정 </h4>

## Introduction
* 주제: 의료 데이터를 이용한 생존 분석 
* 분석 기간: 2023.07 ~ 2023.08
* 데이터: 비공개

<br>

### 1. 폐암 발병 데이터 분석 
<img src="https://github.com/daunJJ/Survival_Analysis/assets/109944763/8abdee83-dcc6-44dc-8bf3-7ea60fb442b3" width="500" height= "250"/>

**생존 상태** 
* 5가지 다중 상태 모형 정의

**적용 모델**
* model1: coxph full model
* model2: 공변량 축소 모델
* model3: 일부 전이가 같다고 가정한 모델

<br>

### 2. 중환자실 호흡기 탈부착 데이터 분석

<img src="https://github.com/daunJJ/Survival_Analysis/assets/109944763/653cbc75-5a1e-4d6f-a6ff-620dd5d11e9e" width="400" height= "200"/>

**생존 상태**
* 순환 전이가 가능한 4가지 다중 상태 모형 정의

**적용 모델**
* 공변량을 사용하지 않은 모델

<br>

<img src="https://github.com/daunJJ/Survival_Analysis/assets/109944763/7cae6567-7305-4ec8-956d-552dc2660734" width="510" height= "230"/>

**생존 상태**
* 순환 전이가 없는 8가지 다중 상태 모형 정의

**적용 모델**
* coxph full model

**추가**
* ventiliation on and off를 공변량으로 한 competing risk 분석
