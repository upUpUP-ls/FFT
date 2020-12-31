using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FFT变换_test
{
    class Program
    {
        /// <summary>
        /// 此处完成的是机器人剧烈磨损，崩刃，以及断刀的判定；过载的话，根据电流的幅值来设置一个经验阈值好了
        /// 在线的怎么弄？
        /// 新建一个队列，设置队列长度为定值，根据采样频率来设，即n值设好了，然后m值根据n值大小也给它写死
        /// 显示的话，显示队列的值可，显示所有的值也可
        /// </summary>
        public static int n = 10000;          //从文件中截取数据的长度,注意n>m
        public static int m = 8192;          //注意必须为2的幂次,与n值最相近的2的幂次
        /// <summary>
        /// 一维频率抽取基2快速傅里叶变换
        /// 频率抽取：输入为自然顺序，输出为码位倒置顺序
        /// 基2：待变换的序列长度必须为2的整数次幂
        /// </summary>
        /// <param name="sourceData">待变换的序列(复数数组)</param>
        /// <param name="countN">序列长度,可以指定[0,sourceData.Length-1]区间内的任意数值</param>
        /// <returns>返回变换后的序列（复数数组）</returns>
        public static Complex[] fft_frequency(Complex[] sourceData, int countN)
        {
            //2的r次幂为N，求出r.r能代表fft算法的迭代次数
            int r = Convert.ToInt32(Math.Log(countN, 2));


            //分别存储蝶形运算过程中左右两列的结果
            Complex[] interVar1 = new Complex[countN];
            Complex[] interVar2 = new Complex[countN];
            interVar1 = (Complex[])sourceData.Clone();

            //w代表旋转因子
            Complex[] w = new Complex[countN / 2];
            //为旋转因子赋值。（在蝶形运算中使用的旋转因子是已经确定的，提前求出以便调用）
            //旋转因子公式 \  /\  /k __
            //              \/  \/N  --  exp(-j*2πk/N)
            //这里还用到了欧拉公式
            for (int i = 0; i < countN / 2; i++)
            {
                double angle = -i * Math.PI * 2 / countN;
                w[i] = new Complex(Math.Cos(angle), Math.Sin(angle));
            }

            //蝶形运算
            for (int i = 0; i < r; i++)
            {
                //i代表当前的迭代次数，r代表总共的迭代次数.
                //i记录着迭代的重要信息.通过i可以算出当前迭代共有几个分组，每个分组的长度

                //interval记录当前有几个组
                // <<是左移操作符，左移一位相当于*2
                //多使用位运算符可以人为提高算法速率^_^
                int interval = 1 << i;

                //halfN记录当前循环每个组的长度N
                int halfN = 1 << (r - i);

                //循环，依次对每个组进行蝶形运算
                for (int j = 0; j < interval; j++)
                {
                    //j代表第j个组

                    //gap=j*每组长度，代表着当前第j组的首元素的下标索引
                    int gap = j * halfN;

                    //进行蝶形运算
                    for (int k = 0; k < halfN / 2; k++)
                    {
                        interVar2[k + gap] = interVar1[k + gap] + interVar1[k + gap + halfN / 2];
                        interVar2[k + halfN / 2 + gap] = (interVar1[k + gap] - interVar1[k + gap + halfN / 2]) * w[k * interval];
                    }
                }

                //将结果拷贝到输入端，为下次迭代做好准备
                interVar1 = (Complex[])interVar2.Clone();
            }

            //将输出码位倒置
            for (uint j = 0; j < countN; j++)
            {
                //j代表自然顺序的数组元素的下标索引

                //用rev记录j码位倒置后的结果
                uint rev = 0;
                //num作为中间变量
                uint num = j;

                //码位倒置（通过将j的最右端一位最先放入rev右端，然后左移，然后将j的次右端一位放入rev右端，然后左移...）
                //由于2的r次幂=N，所以任何j可由r位二进制数组表示，循环r次即可
                for (int i = 0; i < r; i++)
                {
                    rev <<= 1;
                    rev |= num & 1;
                    num >>= 1;
                }
                interVar2[rev] = interVar1[j];
            }
            return interVar2;

        }

        private double[] ChangeQueueToArray(Queue<double> q)
        {
            double[] a = new double[q.Count];
            for (int i = 0; i < q.Count; i++)
            {
                a[i] = Convert.ToDouble(q.Dequeue());
            }
            return a;
        }
        public static void Main(string[] args)
        {
            Complex[] f = new Complex[n];      //存储原始数据
            Complex[] df = new Complex[n];     //存储原始数据的快速傅里叶变换结果
            double[] aa = new double[n];       //存储fft结果的绝对值
            double[] bb = new double[n];       //存储频率区间
            double[] bb1 = new double[n/2];


            //Queue<double> q1 = new Queue<double>(n);
            //Complex[] qq1 = new Complex[n];


            double ab;
            StreamReader sr =File.OpenText("force1.txt");
            string line;
            int p = 0;
            while ((line = sr.ReadLine()) != null)
            {
                string[] value = line.Split(' ');
                foreach (string s in value)
                {
                    if (!string.IsNullOrEmpty(s))
                    {
                        double v = double.Parse(s);
                        f[p] = v;
                        //q1.Enqueue(v);
                        p++;

                        if (p >= n)
                        {
                            break;
                        }
                        if (p > 0 && p < n)
                        {
                            for (int a = p; a < n; a++)
                            {
                                f[a] = 0;
                            }
                        }
                    }
                }
            }

            //for(int i=0;i<n;i++)
            //{
            //    qq1[i] = q1.Dequeue();
            //}
            //此处对复数矩阵f做快速傅里叶变换，结果存放在复数矩阵df中
            df = fft_frequency(f, m);
            //此处对复数矩阵df进行求绝对值计算
            for (int i = 0; i < m; i++)
            {
                Complex c = new Complex(df[i]);
                ab = c.GetModul();
                aa[i] = ab;
            }
            //此处用来计算频率区间
            for(int i=0;i<m;i++)
            {
                bb[i] = i * 10000 / n;
            }
            //bb1[0] = aa[0] / 2048;
            bb1[0] = 0;
            for(int i=1;i<m/2;i++)
            {
                bb1[i] = aa[i] * 2 / n;
            }
            Dictionary<double, double> d = new Dictionary<double, double>();
            //此处将幅值与频率对应起来
            for (int i = 0; i < m/2; i++)
            {
                d.Add(bb1[i], bb[i]);
            }
            List<double> l = new List<double>(d.Keys);
            l.Sort();

            bool t = d.TryGetValue(l[m/2-1], out double row);

            if(l[m / 2 - 1] < 0.005)
            {
                Console.WriteLine("机器人空转");
            }
            else if(row>=166*2-30&&row<=166*2+30)
            {
                Console.WriteLine("机器人正常切削");
            }
            Console.WriteLine("最大幅值对应的频率为："+ row);
            Console.ReadKey();
        }
    }
}
