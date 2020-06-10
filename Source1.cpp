// Bioinformatics - Motif Finding
// by: Anh Dang
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

// print matrix, consensus, and score function
void printSequence(char char_seq[100][100], int num_of_seq, int seq_len)
{
	int acgt[100][4] = { 0 };
	// count ACGT 
	for (int x = 0; x < num_of_seq; x++)
	{
		for (int y = 0; y < seq_len; y++)
		{
			if (char_seq[x][y] == 'A' || char_seq[x][y] == 'a')
			{
				acgt[y][0]++;
			}
			if (char_seq[x][y] == 'C' || char_seq[x][y] == 'c')
			{
				acgt[y][1]++;
			}
			if (char_seq[x][y] == 'G' || char_seq[x][y] == 'g')
			{
				acgt[y][2]++;
			}
			if (char_seq[x][y] == 'T' || char_seq[x][y] == 't')
			{
				acgt[y][3]++;
			}
		}
	}
	// print matrix
	cout << "============" << endl;
	cout << "Input matrix" << endl;
	cout << "============" << endl;
	cout << "A ";
	for (int x = 0; x < seq_len; x++)
	{
		cout << acgt[x][0] << " ";
	}
	cout << endl;
	cout << "C ";
	for (int x = 0; x < seq_len; x++)
	{
		cout << acgt[x][1] << " ";
	}
	cout << endl;
	cout << "G ";
	for (int x = 0; x < seq_len; x++)
	{
		cout << acgt[x][2] << " ";
	}
	cout << endl;
	cout << "T ";
	for (int x = 0; x < seq_len; x++)
	{
		cout << acgt[x][3] << " ";
	}
	cout << endl;
	// print consensus
	int score = 0;
	cout << "===============" << endl;
	cout << "Input consensus" << endl;
	cout << "===============" << endl;
	for (int x = 0; x < seq_len; x ++)
	{
		int temp;
		char consensus;
		if (acgt[x][0] > acgt[x][1])
		{
			temp = acgt[x][0];
			consensus = 'A';
		}
		else { temp = acgt[x][1]; consensus = 'C'; };

		if (temp < acgt[x][2])
		{
			temp = acgt[x][2];
			consensus = 'G';
		}

		if (temp < acgt[x][3])
		{
			temp = acgt[x][3];
			consensus = 'T';
		}
		score += temp;

		cout << consensus << " ";
	}
	//print score
	cout << endl;
	cout << "===========" << endl;
	cout << "Input score = " << score << endl;
	cout << "===========" << endl;
}

// compute TotalDistance
void TotalDistance(string v, string DNA_seq[100], int DNA_num_of_seq)
{
	int DNA_seq_len = DNA_seq[1].length();

	char DNA[100][1000];// max sequence length(l) is 1000
	// copy string array to char array
	cout << "===================" << endl;
	cout << "Input DNA sequences" << endl;
	cout << "===================" << endl;
	for (int x = 0; x < DNA_num_of_seq; x++)
	{
		for (int y = 0; y < DNA_seq_len; y++)
		{
			strcpy_s(DNA[x], DNA_seq[x].c_str());
			cout << DNA[x][y];
		}
		cout << endl;
	}
	char char_v[100];
	
	cout << "User input k-mer: " << endl;

	for (int x = 0; x < v.length(); x++)// copy string to char array
	{
		strcpy_s(char_v, v.c_str());
		cout << char_v[x];
	}
	cout << endl;
	
	// check for similarity between DNA and k-mer
	cout << "===============================================" << endl;
	cout << "Best Match | " << "POS | " << "Hamming Distance" << endl;
	cout << "===============================================" << endl;
	int totalDistance = 0;
	for (int k = 1; k < DNA_num_of_seq; k++)
	{
		int match = 0;
		int temp = 0;
		int bestMatchPos;
		for (int i = 0; i < DNA_seq_len; i++)
		{
			for (int x = 0; x < v.length(); x++)
			{
				if (DNA[k][x + i] == char_v[x])
				{
					match++;
				}
			}
			if (match > temp)
			{
				temp = match;
				bestMatchPos = i;
			}
			match = 0;
		}
		
		for (int x = 0; x < v.length(); x++)
		{
			cout << DNA[k][x + bestMatchPos];
		}
		cout << " | " << bestMatchPos;
		int hammingDist = v.length() - temp;
		cout << "  | " << hammingDist << endl;
		k++;
		totalDistance += hammingDist;
	}
	cout << "================" << endl;
	cout << "Total distance: " << totalDistance << endl;
	cout << "================" << endl;
}

// Brute Force Median String Search 
void bruteForceMedianStr(int k_mer_len, string DNA_seq[100], int DNA_num_of_seq)
{
	int DNA_seq_len = DNA_seq[1].length();

	char DNA[100][1000];// max sequence length(l) is 1000
	// copy string array to char array
	cout << "===================" << endl;
	cout << "Input DNA sequences" << endl;
	cout << "===================" << endl;
	for (int x = 0; x < DNA_num_of_seq; x++)
	{
		for (int y = 0; y < DNA_seq_len; y++)
		{
			strcpy_s(DNA[x], DNA_seq[x].c_str());
			cout << DNA[x][y];
		}
		cout << endl;
	}
	cout << "User input k-mer length: " << k_mer_len << endl;

	int counter[100] = { 0 };
	int temp = 0;
	int pos;
	// first loop
	for (int j = 0; j < DNA_seq_len; j++)
	{
		for (int y = 3; y < DNA_num_of_seq; y++)
		{
			for (int i = 0; i < DNA_seq_len; i++)
			{
				int match = 0;
				for (int x = 0; x < k_mer_len; x++)
				{
					if (DNA[1][x+j] == DNA[y][x + i])
					{
						match++;
					}
				}
				if (match == k_mer_len)
				{
					counter[j]++;
					i = DNA_seq_len;
				}
			}
			y++;
		}

		if (counter[j] > temp)
		{
			temp = counter[j];
			pos = j;
		}
	}

	// second loop
	int counter2[100] = { 0 };
	int temp2 = 0;
	int pos2, pos3;
	int temp3 = 0;
	
	for (int j = 0; j < DNA_seq_len; j++)
	{
		for (int y = 5; y < DNA_num_of_seq; y++)
		{
			for (int i = 0; i < DNA_seq_len; i++)
			{
				int match = 0;
				for (int x = 0; x < k_mer_len; x++)
				{
					if (DNA[3][x + j] == DNA[y][x + i])
					{
						match++;
					}
				}
				if (match == k_mer_len)
				{
					counter2[j]++;
					i = DNA_seq_len;
				}
			}
			y++;
		}

		if (counter2[j] > temp2)
		{
			temp2 = counter2[j];
			pos2 = j;
		}
		if (counter2[j] == temp2)
		{
			temp3 = counter2[j];
			pos3 = j;
		}
	}
	// print motif
	cout << "===================================" << endl;
	cout << "The motif median string of length " << k_mer_len << endl;
	cout << "===================================" << endl;
	if (temp2 < temp)
	{
		for (int x = pos; x < pos + k_mer_len; x++)
		{
			cout << DNA[1][x];
		}
		cout << endl;
	}
	if (temp2 > temp)
	{
		for (int x = pos2; x < pos2 + k_mer_len; x++)
		{
			cout << DNA[3][x];
		}
		cout << endl;
	}
	if (temp3 >= temp)
	{
		for (int x = pos3; x < pos3 + k_mer_len; x++)
		{
			cout << DNA[3][x];
		}
		cout << endl;
	}
}

int main()
{
	// input align txt file
	fstream align;
	align.open("align1.txt");// input file

	string sequence[100];// max number of sequence(n) is 100

	int num_of_seq = 0;
	// copy sequence in string array
	for (int x = 0; x < 100; x++)
	{
		align >> sequence[x];
		if (sequence[x] != "")
		{
			num_of_seq++;
		}
	}
	int seq_len = sequence[0].length();

	char char_seq[100][100];// max sequence length(l) is 100
	// copy string array to char array
	cout << "===============" << endl;
	cout << "Input sequences" << endl;
	cout << "===============" << endl;
	for (int x = 0; x < num_of_seq; x++)
	{
		for (int y = 0; y < seq_len; y++)
		{
			strcpy_s(char_seq[x], sequence[x].c_str());
			cout << char_seq[x][y];
		}
		cout << endl;
	}

	printSequence(char_seq, num_of_seq, seq_len);
	// Input DNA sequences .fa format
	fstream DNA_in;
	DNA_in.open("Sequences.fa.txt");//input DNA sequence .fa file

	string DNA_seq[100];// max number of sequence(n) is 100
	int DNA_num_of_seq = 0;
	// copy sequence in string array
	for (int x = 0; x < 100; x++)
	{
		DNA_in >> DNA_seq[x];
		if (DNA_seq[x] != "")
		{
			DNA_num_of_seq++;
		}
	}
	
	// user input k-mer(v)
	string v;
	cout << "Please input k-mer (Upper-case): ";
	cin >> v;

	TotalDistance(v, DNA_seq, DNA_num_of_seq);

	cout << "Input length of k-mer to find motif median string: ";
	int k_mer_len;
	cin >> k_mer_len;
	bruteForceMedianStr(k_mer_len, DNA_seq, DNA_num_of_seq);

	return 0;
	system("PAUSE");
}