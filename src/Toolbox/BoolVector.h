#pragma once
class BoolVector
{
public:
	BoolVector()
	{
		n = 0;
		vec = nullptr;
	}
	BoolVector(int n_)
	{
		n = n_;
		vec = new bool[n];
	}
	BoolVector(int n_, bool flag_)
	{
		n = n_;
		vec = new bool[n];
		for (int i = 0; i < n; ++i)
			vec[i] = flag_;
	}
	BoolVector(BoolVector& right_)
	{
		if (right_.size() == 0)
		{
			n = 0;
			vec = nullptr;
		}
		else
		{
			n = right_.size();
			vec = new bool[n];
			for (int i = 0; i < n; ++i)
				vec[i] = right_[i];
		}
	}
	BoolVector& operator=(BoolVector &right_)
	{
		if (right_.size() == 0)
			return BoolVector();
		n = right_.size();
		vec = new bool[n];
		for (int i = 0; i < n; ++i)
			vec[i] = right_[i];
		return *this;
	}
	BoolVector& operator=(BoolVector&& right_)
	{
		if (right_.size() == 0)
			return BoolVector();
		n = right_.size();
		vec = right_.vec;
		right_.n = 0;
		right_.vec = nullptr;
		return *this;
	}
	~BoolVector()
	{
		if (vec) { delete vec; vec = nullptr; }
	}
public:
	bool& operator[](int i) { return vec[i]; }
	int size() { return n; }
	void resize(int n_, bool flag_)
	{
		n = n_;
		if (vec)
			delete vec;
		vec = new bool[n];
		for (int i = 0; i < n; ++i)
			vec[i] = flag_;
	}
	void set(bool flag_)
	{
		for (int i = 0; i < n; ++i)
			vec[i] = flag_;
	}
private:
	bool* vec;
	int n;
};