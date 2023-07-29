#include <iostream>
#include <complex>
#include <math.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <unistd.h>
#include <thread>
#include <mutex>
#define pi acos(-1)
using namespace std;


struct Struct
{
	double x, y, z;
	string name;
    int index;
};


Struct hyd[30000];

double box_x, box_y, box_z;

double PBC(double x1, double x2, double y1, double y2, double z1, double z2)
{
	double dert_x = 0;
	double dert_y = 0;
	double dert_z = 0;
	double distance = 0;//����һ������ֵ
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
	distance = (dert_x * dert_x) + (dert_y * dert_y) + (dert_z * dert_z);
	return distance;
}

void PROCESS(int start, int end, int h_n, double x0, double y0, double z0, double dr, std::vector<double>& height_sum, std::vector<double>& count) 
{
	// cout<< "x0 = " << x0 << ", y0 = " << y0 << ", z0 = " << z0 << endl;

	for (int i = start; i < end; i++) //遍历所有该线程内的水合物
	{
		if(hyd[i].z < 9)
		{
			double R = sqrt(PBC(hyd[i].x, x0, hyd[i].y, y0, hyd[i].z, z0));
			int bin_index = int(R/dr); //计算第i个水合物处于哪个区间

			Mutex.lock();
			height_sum[bin_index] += hyd[i].z;
			count[bin_index] ++;
			Mutex.unlock();
		}

	}
}

int main()
{
    int num_threads = 20;
    std::vector<std::thread> threads;
	std::mutex Mutex;
	
	//cout.precision(16);
	
	double a = 0.13458;
	double b = 0.13458;
	string arr;
	ifstream ifs;
	//*************************************************************************************************************
	ifs.open("oneframe.txt", ios::in);
	if (!ifs)
	{
		cout << "*************************************************************************cann't open file" << endl;
		return 0;
	}
	std::vector<string> mol_name;
	std::vector<double> xi;
	std::vector<double> yi;
	std::vector<double> zi;
	std::vector<double> id;
	std::vector<string> type;
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
			string get_name = arr.substr(5, 3);
			string get_type = arr.substr(9, 6);
			string get_id = arr.substr(15, 5);
			string get_x = arr.substr(21, 8);
			string get_y = arr.substr(29, 8);
			string get_z = arr.substr(37, 8);
			mol_name.push_back(get_name);
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

    double x0, y0, z0; //记录质心位置
    ifstream ifs_centro("oneframe_inh_down_centro.txt", ios::in);
    if (ifs_centro.is_open())
    {
        ifs_centro >> x0 >> y0 >> z0;
        ifs_centro.close();
    }
    else
    {
        cout << "Cannot open file" << endl;
    }
    
    // cout<< "x0 = " << x0 << ", y0 = " << y0 << ", z0 = " << z0 << endl;

    int h_n = 0;



	for (int i = 0; i < total; i++)
	{


        if (mol_name[i] == "HYD") //识别所有的HYD分子
		{
			hyd[h_n].x = xi[i];
			hyd[h_n].y = yi[i];
			hyd[h_n].z = zi[i];
            hyd[h_n].index = i;//记录其所处的行数，从0开始，可以直接在vmd里看
			h_n++;
		}

	}

    //对分子的识别没有问题
	//*************************************************************************************************************
    cout<<"h_n = "<<h_n<<endl;
	//*************************************************************************************************************
	
	//清空数组
	vector<double>().swap(xi);
	vector<double>().swap(yi);
	vector<double>().swap(zi);
	vector<double>().swap(id);
	vector<string>().swap(type);

    //#############################################################################################################
    // 创建线程

	// 计算出最大区间数，创造相应大小的数组来存储每个区间的高度总和以及原子数目
	double dr = 0.1;
	int max_bin = int(box_x / dr) + 1;
	cout <<  max_bin << endl;

	vector<double> height_sum(max_bin, 0); //0代表将所有元素初始化为0
	vector<double> count(max_bin, 0);

	int hyd_start, hyd_end;
	cout<< "start\tend\tdert" <<endl;
	cout << "____________________"<<endl;
    for (int i = 0; i < num_threads; i++) //遍历所有线程，对数据进行等分，每份从start-end
	{
        hyd_start = i * (h_n / num_threads);
        hyd_end = (i + 1) * (h_n / num_threads);
		if (i == num_threads - 1) 
		{
        	hyd_end = h_n; //除不尽的时候最后一个线程处理剩余的所有数据
		}
		cout << hyd_start << "\t" << hyd_end << "\t" << hyd_end - hyd_start << endl;
		cout << "____________________"<<endl;
        threads.emplace_back(PROCESS, hyd_start, hyd_end, h_n, x0, y0, z0, dr, height_sum, count);
    }

	//使用一个循环遍历 threads 容器中的所有线程对象，并对每个线程调用 join 方法。这样，主线程会等待所有线程都完成执行后再继续运行。只用一次
    for (auto& t : threads) 
	{
        t.join(); 
    }





	return 0;
}
