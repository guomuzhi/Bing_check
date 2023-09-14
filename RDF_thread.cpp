#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <thread>
#include <mutex>
#define pi acos(-1)
using namespace std;
struct atom
{
	double x,y,z;
	int id;
    int flag;
};
atom wat[50000];
atom hyd[50000];
atom CH4[50000];
double num[1000] = { 0 };
int sum = 0; //记录所有径向区间内原子总数

double box_x, box_y, box_z;


double PBC(double x1, double x2, double y1, double y2, double z1, double z2)
{
	double dert_x = 0;
	double dert_y = 0;
	double dert_z = 0;
	double distance;//����һ������ֵ
	if (abs((x1 - x2)) > (box_x / 2))//x�ں�������
	{
		if (x2 - x1 < 0)
		{
			dert_x = box_x - (x1 - x2);
		}
		else
		{
			dert_x = -(box_x - (x2 - x1));
		}
	}
	else
	{
		dert_x = x2 - x1;
	}
	if (abs((y1 - y2)) > (box_y / 2))//y�ں�������
	{
		if (y2 - y1 < 0)
		{
			dert_y = box_y - (y1 - y2);
		}
		else
		{
			dert_y = -(box_y - (y2 - y1));
		}
	}
	else
	{
		dert_y = y2 - y1;
	}
	if (abs((z1 - z2)) > (box_z / 2))//x�ں�������
	{
		if (z2 - z1 < 0)
		{
			dert_z = box_z - (z1 - z2);
		}
		else
		{
			dert_z = -(box_z - (z2 - z1));
		}
	}
	else
	{
		dert_z = z2 - z1;
	}
	distance = sqrt((dert_x * dert_x) + (dert_y * dert_y) + (dert_z * dert_z));
	return distance;
}


void PROCESS_Judge(int start, int end, int w_n, int h_n)
{
    for (int i = start; i < end; i++)
    {
        if (CH4[i].flag == 0)
        {
            for (int j =0; j < h_n; j++)
            {
                if(PBC(CH4[i].x, hyd[j].x, CH4[i].y, hyd[j].y, CH4[i].z, hyd[j].z) <0.35)
                {
                    CH4[i].flag =1;
                    break;
                }
            }
        }
    }
}


void PROCESS_RDF(int start, int end, int w_n, int h_n, double R, double dr, std::mutex& Mutex)
{
    for (int i = start; i < end; i++)
	{
		for (int j = 0; j < w_n; j++)
		{
			int ceng = 0; //从0层开始
			if (CH4[i].flag == 0)
			{   
				double dis = PBC(CH4[i].x, wat[j].x, CH4[i].y, wat[j].y,  CH4[i].z, wat[j].z);
                Mutex.lock();    
				for (double r = 0; r < R; r += dr)
				{
                    
					if ((dis > r) && (dis <= r + dr))
					{
						num[ceng]++;
						sum++;
					}
					else
					{
						ceng++; //去下一层
					}
				}
                Mutex.unlock();
			}
		}
	}
}


int main()
{
    int num_threads = 10;
    std::vector<std::thread> threads;
    std::mutex Mutex;

	cout.precision(16);
	string arr;
	ifstream ifs;
	ifs.open("panding.gro", ios::in);
	if (!ifs)
	{
		cout << "cann't open file" << endl;
		return 0;
	}
	std::vector<double> xi;
	std::vector<double> yi;
	std::vector<double> zi;
	std::vector<double> id;
	std::vector<string> type;
    std::vector<string> res;
	int total = 0;
	int n_line = 0;

	while (getline(ifs, arr))
	{
		n_line++;
		if (n_line == 2)
		{
			string get_total;
			istringstream totalnumber(arr);
			totalnumber >> get_total;
			total = atof(get_total.c_str());
		}
		if ((n_line > 2) && (n_line <= total + 2))
		{
			string get_res = arr.substr(5, 3);
			string get_type = arr.substr(9, 6);
			string get_id = arr.substr(15, 5);
			string get_x = arr.substr(21, 8);
			string get_y = arr.substr(29, 8);
			string get_z = arr.substr(37, 8);
            res.push_back(get_res);
			id.push_back(atof(get_id.c_str()));
			xi.push_back(atof(get_x.c_str()));
			yi.push_back(atof(get_y.c_str()));
			zi.push_back(atof(get_z.c_str()));
			type.push_back(get_type);
		}
		if (n_line == total + 3)
		{
			string b_x, b_y, b_z;
			istringstream boxsize(arr);
			boxsize >> b_x >> b_y >> b_z;
			box_x = atof(b_x.c_str());
			box_y = atof(b_y.c_str());
			box_z = atof(b_z.c_str());
		}
	}
	ifs.close();
	int w_n = 0;
    int h_n = 0;
    int c_n = 0;

	for (int i = 0; i < total; i++)
	{
		if (type[i] == "   mW ")
		{
			wat[w_n].id = id[i];
			wat[w_n].x = xi[i];
			wat[w_n].y = yi[i];
			wat[w_n].z = zi[i];
			w_n++;
		}
        if (res[i] == "HYD")
		{
			hyd[h_n].id = id[i];
			hyd[h_n].x = xi[i];
			hyd[h_n].y = yi[i];
			hyd[h_n].z = zi[i];
			h_n++;
		}
        if (type[i] == "   mM ")
		{
			CH4[c_n].id = id[i];
			CH4[c_n].x = xi[i];
			CH4[c_n].y = yi[i];
			CH4[c_n].z = zi[i];
            CH4[c_n].flag = 0;
			c_n++;
		}

	}
	// cout << w_n << endl;
    // cout << h_n << endl;
    // cout << c_n << endl;



    // 找到所有液相中甲烷，flag标记为0
    int start, end;
    for (int i = 0; i < num_threads; i++) //遍历所有线程，对数据进行等分，每份从start-end
	{
        start = i * (c_n / num_threads);
        end = (i + 1) * (c_n / num_threads);

		if (i == num_threads - 1) 
		{
        	end = c_n; //除不尽的时候最后一个线程处理剩余的所有数据
		}

        threads.emplace_back(PROCESS_Judge, start, end, w_n, h_n); //传递这两个容器的引用
    }

	//使用一个循环遍历 threads 容器中的所有线程对象，并对每个线程调用 join 方法。这样，主线程会等待所有线程都完成执行后再继续运行。只用一次
    for (auto& t : threads) 
	{
        t.join(); 
    }


    // 对目标组进行RDF计算
    double R = 1;//径向半径nm
	double dr = 0.01;//微分nm

    for (int i = 0; i < num_threads; i++) 
	{
        start = i * (c_n / num_threads);
        end = (i + 1) * (c_n / num_threads);

		if (i == num_threads - 1) 
		{
        	end = c_n; //除不尽的时候最后一个线程处理剩余的所有数据
		}
        cout << start << "\t" << end << "\t" << end - start << endl;
		cout << "____________________"<<endl;

        threads.emplace_back(PROCESS_RDF, start, end, w_n, h_n, R, dr, std::ref(Mutex)); //传递这两个容器的引用
    }



    double density = sum / (4 * R * R * R * pi / 3); //整个径向空间内总数密度
	int ceng = 0;
	ofstream ofs;
	ofs.open("result_RDF.csv");
	for (double r = 0; r < R; r += dr)
	{
		ofs << (r + dr) << "," << num[ceng] / (density * (4 * pi * (r + dr / 2) * (r + dr / 2) * dr)) << endl;//归一化处理，每层的数目/平均数密度*该层的体积
	}
	ofs.close();

	return 0;
}

