## leiyuyu的ETRU说明
**写在前面**
过程中出了很多版本，也许这个版本还不太完善。先这样吧。过些天闲下来，再仔细校对很多问题，有时候可能在修改过程中不小心把有的地方改了。其次可能是没有分开编写，一整个cpp导致的问题，不过咱们尽管先占个坑位。
**关于接口**

代码位于group_7\group_7\group_7.cpp；也可以自己新建控制台项目，直接把文件group_7.cpp的内容拷进去，执行，这样最省事。
**关于接口**

cpp里我写了很多测试函数，只需要把main里面的测试函数换掉，即可运行测试函数；
测试包括：
对密钥生成时间的测试：
`keyGen_only()`
对多项式求逆以及三种乘法的验证：
`test_inverse_and_multiplications()`
对密文加密解密的验证：
`interface_with_data()`
对输入信息进行多项式转化再转回来：
`test_translation() `
爱森斯坦数和虚数之间的转换：
`test_EI_to_Complex()`

**关于成功与否**

如果只是密钥生成的话，N一般可以在800以上不报错，偶尔会越界错误，可以再走一次；

数据加解密的时候，可能调用链长的缘故，有时候会不打印解密还原后的信息。鉴于这个情况，在最后的阶段，我把明文多项式和密文对应的多项式都打印出来了，可以自己试一下。

**关于测试不同乘法**

这里其实保留了三个乘法：

```
multiplication_normal(), multiplication_convolution(), multiplicaiton_FFT()
```

如果想比较运行时间的话，可以在cpp文件中Ctrl+F搜索：
`newer = old.subtraction(quotientAndRemainder.first.multiplication_convolution(neww, p), p);`该行语句落位在这个函数：`polyInverse(Poly P, EI p)`内；

只要把`multiplication_convolution`换成`multiplication_FFT`，或者换成`multiplication_normal，然后在main中运行`keyGen_only()`，即可观察不同算法的求逆时间，需要指出的是normal算法下，N不宜大（<100），否则会运行好长时间出不来.FFT的算法是初步的版本，只有在达到N在600+的时候才能感受到很明显的加速效果。

**关于更多改进**
在IFFT之前先进行截断，让系数变短，是另外一个盛时的方面，我有一点具体的想法了，需要一些时间来尝试；其次我想把静态数组换成vector，之前试的时候试错了，以为不可以做到呢。

