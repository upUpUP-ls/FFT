# FFT
使用C#编写的快速傅里叶变换小demo，处理加工过程中三相力信号和电流信号，变换之后的处理过程其实就是参照matlab软件中fft函数对信号的处理过程
只是我稍微做了一些修改，利用C#的dictionary绑定频率和幅值，可以提取最大幅值对应的频率，或者说可以提取第几大幅值对应的频率，这些频率大小其
实可以对应铣刀加工过程中的齿频，就可以作为判断加工过程中刀具状态的一个标准
