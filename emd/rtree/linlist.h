#ifndef __LINLIST
#define __LINLIST

#include <stdio.h>

// LinList  (SLink)

struct Linkable {
public:
	int son,dimension,level;
	float *bounces;
	float distanz;

	Linkable(int dim) { 
		dimension = dim;
	  	bounces = new float[2 * dim];
    }

	~Linkable() { 
		delete [] bounces;
	}
};

struct SLink {
    Linkable *d;
    SLink *next,*prev;
    SLink();
    ~SLink();
};

// LinList

class LinList {
protected:
    SLink *first;         // Rootzeiger des Datenbestands
    SLink *last;          // Zeiger auf letztes Element
    int anz;                    // Anzahl der belegten Elemente in der Liste
    SLink *akt;           // zeigt auf aktuelles Element
    
public:
	int akt_index;              // Index des zuletzt mit get geholten Elements
    LinList();
    virtual ~LinList();
    int get_num()               // gibt Anzahl der im Index belegten Elements
        { return anz; }         // zurueck
	
    void check();               // ueberprueft Konsistenz der Liste
    void print();

    void insert(Linkable *f);       // haengt ein Element vorne an die Liste an
    bool erase();               // loescht aktuelles Element aus der Liste

    Linkable * get(int i);          // liefert i-tes Element
    Linkable * get_first();         // liefert erstes Element im Index
    Linkable * get_last();          // liefert erstes Element im Index
    Linkable * get_next();          // liefert naechstes Element im Index
    Linkable * get_prev();          // liefert vorhergehendes Element im Index
};

#endif  // __LINLIST
