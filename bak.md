# 待优化 
* 优化产生分布的代码，主要是优化快速寻找reads上微卫星边界的方法
* 产生的分布的异常的排除，特别是异常比对，或者segmental duplication的情况
    * 通过1.5*IQR来过滤 
      * IQR=(Q3-Q1) 如果IQR<20 IQR=20
      * 去除<Q1-QR*1.5 或者 >Q3+IQR*1.5
