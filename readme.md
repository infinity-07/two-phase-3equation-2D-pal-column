# 有限体积法求解双曲型方程

## 已经做到了的
- 使用 WENO-JS 格式实现五阶精度（对流，欧拉，两相流）
- 周期边界条件
- 实现了两相流的保界限制器

## 仍然需要优化的：
- 增加两种重构方式，
- Zhang, R., Zhang, M., & Shu, C.-W. (2015). On the Order of Accuracy and Numerical Performance of Two Classes of Finite Volume WENO Schemes. Communications in Computational Physics, 9(3), 807-827. https://doi.org/10.4208/cicp.291109.080410s 

- 增加 thinc
- 研究如何输出成 paraview 格式
- 增加断点保存功能
- 增加设置是否输出虚拟单元的config
- WENO 重构调整为输出一个数组

- lambda delta x 超过上界的时候，增加一个thinc呢，增加一个 lambda delta x 指示器

GarnetMccarvilleuuUcV@gmail.com --- acwhqd9TeoBI6 --- GarnetMccarvilleuuUcV@outlook.com

选择性的使用保界限制器，如果theta=1就不应用


增加一个中间时间节点输出

增加一个每次输出数据结果的时候文件夹名字修改为 testcase+numX+numY+xxxxxxxx

搞清楚为什么修改了 getPressure 之后速度会变

检查边界处理情况

增加 输出各线程的结果
增加 输出ghost单元的结果



仔细检查一遍导致 UPV 不稳定的原因

1. 除以一个小数
2. 黎曼求解器
3. 为什么更换了压力的求解结果速度更稳定了
4. 为什么 eps 调小就稳定

重新增加一个变量，边界的单元平均值，方便和其他的重构方法配合


增加了无量纲化

只有z1使用指数

最好还能启动一次，完成所有的算例

设置testcase 的时候把z1改小

存在的问题：
HLLC有震荡
PP-limiter位置不对，
边界有震荡


检查边界条件设置是否正确，检查非线性权重设置方案