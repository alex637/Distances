#include <algorithm>


template <typename T, typename Comp>
class Heap
{
public:
	int size, n;
	T * arr;
	
	
	Heap(int capacity, int * indexes = nullptr)
	{
		size = capacity;
		n = 0;
		indices = indexes;
		T * x = new T[size];
		arr = x;
	}

private:
	int * indices;
	Comp comp;
	
	void SiftDown(int i)
	{
		int left = 2 * i + 1, right = left + 1;
		int index_min = i;
		if (left < n && comp(arr[i], arr[left]))
			index_min = left;
		if (right < n && comp(arr[i], arr[right]))
			index_min = right;
		if (index_min != i)
		{
			std::swap(arr[i], arr[index_min]);
			if (indices != nullptr)
				std::swap(indices[arr[i].vertex], indices[arr[index_min].vertex]);
			SiftDown(index_min);
		}
	}

	void SiftUp(int i)
	{
		while (i > 0)
		{
			int parent = (i - 1) / 2;
			if (comp(arr[i], arr[parent]) || arr[i] == arr[parent])
				return;
			std::swap(arr[i], arr[parent]);
			if (indices != nullptr)
				std::swap(indices[arr[i].vertex], indices[arr[parent].vertex]);
			i = parent;
		}
	}

public:
	~Heap()
	{
		delete[] arr;
	}

	bool empty() { return n == 0; }

	void push(T x)
	{
		
		n++;
		if (n > size)
			return;
		arr[n - 1] = x;
		if (indices != nullptr)
			indices[x.vertex] = n - 1;
		SiftUp(n - 1);
	}

	T pop()
	{
		n--;
		//if (n < 0) return;
		T tmp = arr[0];
		arr[0] = arr[n];
		SiftDown(0);
		if (indices != nullptr)
			indices[tmp.vertex] = -1; // the element has been deleted
		return tmp;
	}

	void refresh(int index)
	{
		SiftUp(index);
		SiftDown(index);
	}
};