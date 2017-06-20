#include<iostream>
#include<fstream>  
#include<string>
#include<vector>
#include<time.h>
#include<cmath>
#include<algorithm>


using namespace std;
//Function to find the rows of file
int get_rows(string fila_name)
{
	int n_i = 0;
	ifstream infile(fila_name);//此处默认的文件打开方式为“以输出的方式打开”
	char str[10000];//N是定义的常数，目的是为了读取足够长的行

	vector<double> prices;
	while (!infile.eof())
	{
		infile.getline(str, sizeof(str));//此处默认的终止标识符为‘\n’
		n_i++;
	}
	return n_i - 1;
}

//Function to read txt into 2-D vector
vector< vector<string> > Read_2_string(int rows, int cols, string file_name)
{
	vector<vector<string>> data(rows, vector<string>(cols));
	ifstream myReadFile;
	myReadFile.open(file_name);

	while (!myReadFile.eof()) {

		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				myReadFile >> data[i][j];
			}
		}

	}
	return data;
}

//Function to copy certain rows of a matrix
vector< vector<double> > Copy_matrix(vector< vector<double> >data, int start_ps, int end_ps, int num_stk)
{
	//Select the matrix
	int len = end_ps - start_ps + 1;
	vector< vector<double> > copy_matrix(len, vector<double>(num_stk));

	for (int i = start_ps; i <= end_ps; i++)
	{
		copy_matrix[i - start_ps] = data[i];
	}

	return copy_matrix;

}


//Function to estimate the beta given time period
vector<double> Get_beta(vector< vector<double> >stk_return, vector<double> mkt_return, \
						vector<double> rf, int start_ps, int end_ps)
{
	int num_stk = stk_return[0].size();
	
	int n = end_ps - start_ps + 1;
	vector<double> beta(num_stk);
	//Stock by stock
	for (int s = 0; s < num_stk; s++)
	{
		//Initialize the sums at the beginning of every stock 
		double x, y, sum_x = 0, sum_y = 0, sum_x_2 = 0, sum_xy = 0;
		//Day by day
		for (int i = start_ps; i <= end_ps; i++)
		{
			x = mkt_return[i] - rf[i];
			y = stk_return[i][s] - rf[i];
			sum_x += x;
			sum_y += y;
			sum_x_2 += x * x;
			sum_xy += x * y;
		}
		//OLS estimation, calculate the coefficient
		beta[s] = (n * sum_xy - sum_x * sum_y) / (n * sum_x_2 - sum_x * sum_x);
	}
	return beta;
}

//Get the rolling beta
vector< vector<double> > Get_rolling_beta(int roll_window, vector< vector<double> >stk_return, \
											vector<double> mkt_return, vector<double> rf)
{
	int time_periods = mkt_return.size() - roll_window;
	int num_stk = stk_return[0].size();
	vector< vector<double> > rolling_beta(time_periods, vector<double>(num_stk));
	for (int i = 0; i < time_periods; i++)
	{
		int start_ps = i, end_ps = i + roll_window - 1;
		rolling_beta[i] = Get_beta(stk_return, mkt_return, rf, start_ps, end_ps);
	}
	return rolling_beta;
}

//Function to estimate the Kalman beta
vector< vector<double> > Get_Kal_beta(vector< vector<double> >stk_return, vector<double> mkt_return, \
							vector<double> rf, double Q, double R, double P0, vector<double> init_beta)
{
	int num_days = stk_return.size();
	int num_stk = stk_return[0].size();
	vector< vector<double> > Kal_beta(num_days, vector<double>(num_stk));
	double A = 1.0, K0 = 1;
	double beta_minus, beta_hat, P_minus, P_hat, K, epslon;
	double beta_0;
	for (int s = 0; s < num_stk; s++)
	{
		beta_0 = init_beta[s];
		beta_minus = A * beta_0;
		P_minus = A * A * P0 + Q;
		beta_hat = beta_0;
		P_hat = P0;
		K = K0;
		Kal_beta[0][s] = beta_hat;
		for (int d = 1; d < num_days; d++)
		{
			//Calculate the paremeters
			beta_minus = A * beta_hat;
			P_minus = pow(A, 2) * P_hat + Q;
			epslon = stk_return[d][s] - rf[d] - beta_minus * (mkt_return[d] - rf[d]);
			//Update the parameters
			K = (P_minus * (mkt_return[d] - rf[d])) / \
				(P_minus * pow((mkt_return[d] - rf[d]), 2) + R);
			beta_hat = beta_minus + K * epslon;
			P_hat = pow((1 - K * (mkt_return[d] - rf[d])), 2) * P_minus + pow(K, 2) * R;
			//Store the corrected beta estimation
			Kal_beta[d][s] = beta_hat;
		}
	}
	return Kal_beta;
}

//Function to print the data
void Print_data(vector< vector<string> > data)
{
	int rows = data.size();
	int cols = data[0].size();
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			cout << data[i][j] << " ";
		}
		cout << endl;
	}
}

void Print_data(vector< vector<double> > data)
{
	int rows = data.size();
	int cols = data[0].size();
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			cout << data[i][j] << " ";
		}
		cout << endl;
	}
}

//Function to print the vector
void Print_vector(vector<string> data)
{
	for (int i = 0; i < data.size(); i++)
	{
		cout << data[i] << endl;
	}
}

void Print_vector(vector<double> data)
{
	for (int i = 0; i < data.size(); i++)
	{
		cout << data[i] << endl;
	}
}

//Function to test the Kalman beta strategy
vector< vector<double> > Back_test_beta(vector< vector<double> >stk_return, vector<double> mkt_return, \
	vector<double> rf, double Q, double R, double P0, vector<double> init_beta)
{
	int num_days = stk_return.size();
	int num_stk = stk_return[0].size();
	vector< vector<double> > Kal_beta(num_days, vector<double>(num_stk));
	double A = 1.0, K0 = 1;
	double beta_minus, beta_hat, P_minus, P_hat, K, epslon;
	double beta_0;
	for (int s = 0; s < num_stk; s++)
	{
		beta_0 = init_beta[s];
		beta_minus = A * beta_0;
		P_minus = A * A * P0 + Q;
		beta_hat = beta_0;
		P_hat = P0;
		K = K0;
		Kal_beta[0][s] = beta_hat;
		for (int d = 1; d < num_days; d++)
		{
			//Calculate the paremeters
			beta_minus = A * beta_hat;
			P_minus = pow(A, 2) * P_hat + Q;
			epslon = stk_return[d][s] - rf[d] - beta_minus * (mkt_return[d] - rf[d]);
			//Update the parameters
			K = (P_minus * (mkt_return[d] - rf[d])) / \
				(P_minus * pow((mkt_return[d] - rf[d]), 2) + R);
			beta_hat = beta_minus + K * epslon;
			P_hat = pow((1 - K * (mkt_return[d] - rf[d])), 2) * P_minus + pow(K, 2) * R;
			//Store the corrected beta estimation
			Kal_beta[d][s] = beta_hat;
		}
	}
	return Kal_beta;
}

//Function to write into csv file
void Write_to_csv(string csv_file, vector< vector<double> > data)
{
	ofstream myfile;
	myfile.open(csv_file);
	int rows = data.size();
	int cols = data[0].size();
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (j < cols - 1)
			{
				myfile << data[i][j] << ",";
			}
			else
			{
				myfile << data[i][j] << "\n";
			}
		}

	}
	myfile.close();
}

int main()
{
		//Running starts
		clock_t start, finish;
		double totaltime;
		start = clock();
	
		//Tell the file name
		string file_prices = "Last_price_small.txt";
		const int rows = 1494, cols = 5;
	
		//Read txt into string matrix
		vector< vector<string> > raw_data = Read_2_string(rows, cols, file_prices);
		
	
	
		//String vector to store the days
		vector<string> days;
		for (int i = 0; i < raw_data.size(); i++)
		{
			days.push_back(raw_data[i][0]);
		}
	
		//Establish a matrix to store the stock returns
		vector< vector<double> > stk_return(rows - 1, vector<double>(cols - 1));
		for (int i = 0; i < rows - 1; i++)
		{
			for (int j = 0; j < cols - 1; j++)
			{
				//Calculate the stock daily return
				stk_return[i][j] = stod(raw_data[i + 1][j + 1]) / stod(raw_data[i][j + 1]) - 1;
			}
		}
		int num_stk = stk_return[0].size();
		int num_days = rows - 1;
	
		//Read in market price and risk-free return
		string file_market = "Market_rf_1.txt";
		vector< vector<string> > market_rf = Read_2_string(rows, 3, file_market);
	
		//Vector to store the market return
		vector<double> mkt_return;
		//Vector to store the risk-free rate
		vector<double> rf;
		for (int i = 0; i < num_days; i++)
		{
			mkt_return.push_back(stod(market_rf[i + 1][1]) / stod(market_rf[i][1]) - 1);
			rf.push_back(stod(market_rf[i + 1][2]));
		}
	
		//Select the first n_test days to do the prior estimating of parameters
		int n_test = 100;
	
		//Calculate initial beta for each stock
		vector<double> init_beta(num_stk);
		init_beta = Get_beta(stk_return, mkt_return, rf, 0, n_test);
		//cout << "Initially the betas were:" << endl;
	
		////Get the rolling beta
		//int roll_window = 100;
		//vector< vector<double> > rolling_beta = Get_rolling_beta(roll_window, stk_return, mkt_return, rf);
		// 
		////Print the rolling beta
		//cout << "Rolling beta is of size [" << rolling_beta.size() << "," \
			//		<< rolling_beta[0].size() << "] as follows: \n" << endl;
////Print_data(rolling_beta);
//cout << "***************************************************************************************" <<endl;

//Kalman estimate betas
//First input the parameters for correction
	double Q = 0.002, R = 0.003, P0 = 0.0025;
	//Then input the data range outside the initial ones
	int rest_days = num_days - n_test;
	vector<double> kal_mkt(rest_days), kal_rf(rest_days);
	//Select the two vectors
	copy(mkt_return.begin() + n_test, mkt_return.end(), kal_mkt.begin());
	copy(rf.begin() + n_test, rf.end(), kal_rf.begin());
	//Select the matrix
	vector< vector<double> > kal_stk_return = Copy_matrix(stk_return, n_test, num_days - 1, num_stk);

	//Kalman estimation of betas
	vector< vector<double> > Kal_beta = Get_Kal_beta(kal_stk_return, kal_mkt, kal_rf, Q, R, P0, init_beta);
	//Print the Kalman beta
	/*cout << "Kalman beta is of size [" << Kal_beta.size() << "," \
	<< Kal_beta[0].size() << "] as follows: \n" << endl;
	Print_data(Kal_beta);*/


	//Test beta
	//Select part of the different datasets
	int begg = 1, endd = 20;

	/*vector< vector<double> > small_beta = Copy_matrix(Kal_beta, begg, endd, num_stk);
	cout << "Small sample of beta with size:" << small_beta.size() << ", " << small_beta[0].size() << endl;
	Print_data(small_beta);
	cout << endl;*/

	Print_data(stk_return);

	//Select part of the stock return data
	vector< vector<double> > small_stk = Copy_matrix(stk_return, begg, endd, num_stk);
	//cout << "\nSmall sample of stock return size:" << small_stk.size() << ", " << small_stk[0].size() << endl;
	//Print_data(small_stk);
	//cout << endl;
	//cout << endl;

	//Print_vector(mkt_return);
	//Select part of market return data
	int s_days = endd - begg + 1;
	vector<double> small_mkt(s_days);
	copy(mkt_return.begin() + begg, mkt_return.begin() + endd + 1, small_mkt.begin());
	//cout << "\nSmall sample of market return size:" << small_mkt.size() << endl;
	//Print_vector(small_mkt);

	//Find beta signals
	vector<double> port_return(s_days);




	////Write to csv-Rolling beta and Kalman beta
	//string Rolling_beta_csv = "Rolling_beta.csv";
	//string Kalman_beta_csv = "Kalman_beta.csv";
	//Write_to_csv(Rolling_beta_csv, rolling_beta);
	//Write_to_csv(Kalman_beta_csv, Kal_beta);

	finish = clock();
	totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "\nTime used is " << totaltime << "s" << endl;


	return 0;
}