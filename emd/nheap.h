/* used by heap.c */

#ifndef NHEAP_HH
#define NHEAP_HH

#define MAXELEMS 1000

using namespace std;

template <class V>
class fHeap {
public:
    typedef struct el
    {
        V value;
        int idx;
    } elem;

    struct elemcomp {
        bool operator () (elem left, elem right) const
        { return left.value<right.value; }
    };

    enum sign {INC=0, DEC=1};

    /*returns the parent of a heap position*/
//    inline int parent(int posel) {
//        if (posel%2) return posel/2;
//        else return (posel-1)/2;
//    };

    inline int parent(int posel) {
    	return (posel-1) >> 1;
    }

    fHeap(int sign=INC) {
		num_elems = 0;
        this->sign=sign;
    };

    ~fHeap() {
        vector<elem>().swap(heap);
        vector<int>().swap(order);
    };

	void clear() {
		heap.clear();
		order.clear();
		num_elems = 0;
	}

	bool isExisted(int idx, V &v) {
		if (num_elems>0 && idx<order.size() && order[idx]!=-1) {
			v=heap[order[idx]].value;
			return true;
		} else
			return false;
	};

	bool isExisted(int idx) {
		if (num_elems>0 && idx<order.size() && order[idx]!=-1) {
			return true;
		} else
			return false;
	};

	// old remove function with bug!!!

	//bool remove(int idx) {
	//	if (num_elems>0 && idx<order.size() && order[idx]!=-1) {
	//		int posel = order[idx];
 //           heap[posel].value = -1;
 //           moveup(posel);

 //           V value;
 //           dequeue(idx, value);

	//		return true;
	//	}
	//	return false;
	//};


	bool remove(int idx)
	{
		if (num_elems > 0 && idx < order.size() && order[idx] != -1)
		{
			int pos = order[idx];
			order[idx] = -1;
			if (num_elems - 1 > 0)
			{
				order[heap[num_elems-1].idx] = pos;
				heap[pos].idx = heap[num_elems-1].idx;
			//heap[pos].value = heap[endpos].value;
				num_elems--;
				updatequeue(heap[pos].idx, heap[num_elems].value);
			}
			else num_elems--;
			return true;
		}
		return false;
	}

    void enqueue(int idx, V value) {
		//if (num_elems>0 && idx<order.size() && order[idx]!=-1) {
        //    updatequeue(idx, value);
        //    return;
        //}

        V tmp;
        int tmpidx,p,posel;

        posel = num_elems; //last position
        if (heap.size()<num_elems+1)
            heap.push_back(elem());
        heap[num_elems].value = value;
        heap[num_elems++].idx = idx;
        if (idx>=order.size())
            order.insert(order.end(), idx-order.size()+1, -1);
        order[idx] = posel;

        while (posel >0) {
            p=parent(posel);
            //if (value<heap[p].value) {
            if (compare(value, heap[p].value)) {
                /* swap element with its parent */
                tmp = heap[p].value;
                heap[p].value = value;
                heap[posel].value = tmp;

                tmpidx = heap[p].idx;
                heap[p].idx = idx;
                heap[posel].idx = tmpidx;

                /* adjust order */
                order[idx] = p;
                order[tmpidx] = posel;

                posel = parent(posel);
            } else
                break;
        }
    };

	int dequeue(int &idx) {
		if (num_elems==0) /* empty queue */
            return 0;

        //value = heap[0].value;
		idx = heap[0].idx;
		//printf("dequeue %d", idx);
        order[idx] = -1; /* not in queue */
		if (num_elems-1>0) {
			heap[0].value = heap[num_elems-1].value;
			heap[0].idx = heap[num_elems-1].idx;
			order[heap[0].idx] = 0;
			num_elems--;
			movedown(0);
		} else
			num_elems--;

        return 1;
    };

	int dequeue() {
        if (num_elems==0) /* empty queue */
            return 0;

        //value = heap[0].value;
		int idx = heap[0].idx;
		//printf("dequeue %d ", idx);
        order[idx] = -1; /* not in queue */
		if (num_elems-1>0) {
			heap[0].value = heap[num_elems-1].value;
			heap[0].idx = heap[num_elems-1].idx;
			order[heap[0].idx] = 0;
			num_elems--;
			movedown(0);
		} else
			num_elems--;

        return 1;
    };

    int dequeue(int &idx, V &value) {
        if (num_elems==0) /* empty queue */
            return 0;

        value = heap[0].value;
        idx = heap[0].idx;
        order[idx] = -1; /* not in queue */
		if (num_elems-1>0) {
			heap[0].value = heap[num_elems-1].value;
			heap[0].idx = heap[num_elems-1].idx;
			order[heap[0].idx] = 0;
			num_elems--;
			movedown(0);
		} else
			num_elems--;

        return 1;
    };

    void updatequeue(int idx, V new_val) {
		int posel;

        posel = order[idx];

        if (new_val==heap[posel].value)
            return;

        //if (new_val<=heap[posel].value)
        if (!compare(heap[posel].value, new_val))
        {
            heap[posel].value = new_val;
            moveup(posel);
        } else
        {
            heap[posel].value = new_val;
            movedown(posel);
        }
    };

    void movedown(int md_start) {
        V tmp;
        int tmpidx;
        int posel = md_start; //root
        int swap;
        /*while posel is not a leaf and heap[posel].value > any of childen*/
        while ((posel<<1)+1 < num_elems) { // exists a left son
            if ((posel<<1)+2 < num_elems) { // exists a right son
                //if (heap[posel*2+1].value<heap[posel*2+2].value) // find minimum son
                if (compare(heap[(posel<<1)+1].value, heap[(posel<<1)+2].value))
                    swap = (posel<<1)+1;
                else
                    swap = (posel<<1)+2;
            } else
                swap = (posel<<1)+1;

            //if (heap[posel].value > heap[swap].value) { /*larger than smallest son*/
            if (compare(heap[swap].value, heap[posel].value)) {
                /*swap elements*/
                tmp = heap[swap].value;
                heap[swap].value = heap[posel].value;
                heap[posel].value = tmp;

                tmpidx = heap[swap].idx;
                heap[swap].idx = heap[posel].idx;
                heap[posel].idx = tmpidx;

                /*swap orders*/
                order[heap[swap].idx] = swap;
                order[heap[posel].idx] = posel;

                posel = swap;
            } else {
                break;
            }
        }
    };

    void moveup(int md_start=0) {
        V tmp;
        int tmpidx;
        int posel = md_start; //root
        int p;
        int idx=heap[posel].idx;

        /*while posel is not a root and heap[posel].value < parent */
        while (posel >0) {
            p=parent(posel);
            //if (heap[posel].value<heap[p].value) {
            if (compare(heap[posel].value,heap[p].value)) {
                /* swap element with its parent */
                tmp = heap[p].value;
                heap[p].value = heap[posel].value;
                heap[posel].value = tmp;

                tmpidx = heap[p].idx;
                heap[p].idx = heap[posel].idx;
                heap[posel].idx = tmpidx;

                /* adjust order */
                order[heap[p].idx] = p;
                order[heap[posel].idx] = posel;

                posel = p;
            } else
                break;
        }
    };

    void getTop(int &idx, V &value) {
        idx=heap[0].idx;
        value=heap[0].value;
    };

    void reset() { 
        heap.clear();
        order.clear();
        num_elems=0;
    };

    void reset2() { 
        heap.clear();
        num_elems=0;
    };

    V getTopValue() { return heap[0].value; };
	int getTopIdx() { return heap[0].idx; };
	V getSecondTopValue() { 
		return heap[1].value < heap[2].value ? heap[1].value : heap[2].value;
	}

	int size() { return num_elems; }

    void print_heap() {
//        for (int i=0; i<num_elems; i++)
//            printf("(%d) ", order[heap[i].idx]);
//        printf("\n");
		//system("pause");
    };

//private:
    bool compare(V a, V b) { return (sign==INC)? (a<b) : (a>b); };

    vector<elem> heap;
    vector<int> order;
    int num_elems;
    int sign;
};

///
/// Heap - end

#endif
