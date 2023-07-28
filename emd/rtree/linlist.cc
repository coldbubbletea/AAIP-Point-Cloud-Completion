#include "linlist.h"
#include "gendef.h"

SLink::SLink() {
    d = NULL;
    next = prev = NULL;
}

SLink::~SLink() {
    delete d;
}

LinList::LinList() {
    anz = 0;
    akt_index = -1;
    akt = first = last = NULL;
}

LinList::~LinList() {
    SLink *hf;
    akt = first;
    while (akt != NULL) {
		hf = akt->next;
		delete akt;
		akt = hf;
    }
}

void LinList::insert(Linkable *f) {
    SLink *sd;

    // neuen Rahmen erzeugen
    sd = new SLink;
    sd->d = f;

    // Zeiger umbiegen
    sd->next = first;
    sd->prev = NULL;
    if (first != NULL)
	first->prev = sd;

    // first, last und anz berichtigen
    anz++;
    first = sd;
    if (last == NULL)
	last = sd;

    // Position ist undefiniert
    akt = NULL;
    akt_index = -1;
}

bool LinList::erase() {
    SLink *n_akt;

    // Liste leer oder akt nicht definiert?
    if (akt)
    {
	// Element ist erstes Element
	if (akt == first)
	{
	    // Element ist einziges Element
	    if (akt == last)
	    {
		akt_index = -1;
		first = last = NULL;
		n_akt = NULL;
	    }
	    else
	    {
		// Element ist erstes Element, aber nicht leztes
		(akt->next)->prev = NULL;
		first = akt->next;
		n_akt = first;
		akt_index = 0;
	    }
	}
	else
	{
	    // Element ist letztes Element
	    if (akt == last)
	    {
		(akt->prev)->next = NULL;
		last = akt->prev;
		n_akt = NULL;
		akt_index = -1;
	    }
	    else
	    // Element ist mitten in der Liste
	    {
		(akt->next)->prev = akt->prev;
		(akt->prev)->next = akt->next;
		n_akt = akt->next;
		akt_index++;
	    }
	}

	// Speicher freigeben
	delete akt;

	// aktuelles Element setzen
	akt = n_akt;

	// anz berichtigen
	anz--;
	return true;
    }
    return false;
}

Linkable* LinList::get(int i) {
    bool ahead;   // wenn ahead true ist, wird in next-Richtung gesucht
    int j;

    // liegt das i-te Element ueberhaupt in der Liste?
    if (i >= anz)
	return NULL;

    // ist die Liste schon auf das i-te Element positioniert?
    if (i == akt_index)
	return akt->d;

    // hat eine Positionierung der Liste stattgefunden?
    if (akt_index == -1)
    {
	// i liegt naeher an first, als an last
	if (i < (anz / 2))
	{
	    akt = first;
	    akt_index = 0;
	    ahead = true;
	}
	else
	{
	    akt = last;
	    akt_index = anz - 1;
	    ahead = false;
	}
    }
    else
    {
	// die gewuenschte Position liegt vor der aktuellen
	if (i < akt_index)
	{
	    // liegt i naeher an first, als an akt_index?
	    if ((akt_index - i) > i)
	    {
		akt = first;
		akt_index = 0;
		ahead = true;
	    }
	    else
		ahead = false;
	}
	else
	{
	    // liegt i naeher an last, als an akt_index?
	    if ((i - akt_index) > ((anz-1) - i))
	    {
		akt = last;
		akt_index = anz - 1;
		ahead = false;
	    }
	    else
		ahead = true;
	}
    }
    
    if (ahead)
    {
	for (j = akt_index; j < i; j++)
	{
	    if (!akt)
		error("LinList::get: List seems to be inkonsistent", true);
	    akt = akt->next;
	}
    }
    else
    {
	for (j = akt_index; j > i; j--)
	{
	    if (!akt)
		error("LinList::get: List seems to be inkonsistent", true);
	    akt = akt->prev;
	}
    }
    akt_index = i;
    return akt->d;
}

Linkable* LinList::get_first()
{
    akt = first;

    if (akt != NULL)
    {
	akt_index = 0;
	return akt->d;
    }
    else
	return NULL;
}

Linkable* LinList::get_last() {
    akt = last;

    if (akt != NULL)
    {
	akt_index = anz - 1;
	return akt->d;
    }
    else
	return NULL;
}

Linkable* LinList::get_next() {
    akt = akt->next;
    if (akt != NULL)
    {
	akt_index++;
	return akt->d;
    }
    else
    {
	akt_index = -1;
	return NULL;
    }
}

Linkable* LinList::get_prev()
{
    akt = akt->prev;

    if (akt != NULL)
    {
	akt_index--;
	return akt->d;
    }
    else
    {
	akt_index = -1;
	return NULL;
    }
}

void LinList::print()
{
    SLink *cur = first;

    while (cur != NULL)
    {
	  printf("%d %f %f %f %f\n", cur->d->son, cur->d->bounces[0], cur->d->bounces[1], cur->d->bounces[2],  cur->d->bounces[3]);
	  cur = cur->next;
    }
}

void LinList::check()
{
    SLink *f, *old_f;
    int myanz;
    char buffer[255];

    old_f = first;
    // Liste muss ganz leer sein
    if (old_f == NULL)
    {
	if (last != NULL)
	    error("LinList::check: first == NULL, last != NULL", false);
	if (anz != 0)
	    error("LinList::check: first == NULL, anz != 0", false);
	return;
    }

    myanz = 1;
    if (old_f->prev != NULL)
    {
	error("LinList::check: Listenkopf.prev ungleich NULL", false);
	return;
    }

    for (f = old_f->next; f != NULL; f = f->next)
    {
	if (f->prev != old_f)
	{
	    error("LinList::check: Rueckwaertsverkettung fehlerhaft", false);
            return;
        }
	if (old_f->next != f)
	{
	    error("LinList::check: Vorwaertsverkettung fehlerhaft", false);
            return;
        }
	old_f = f;

	myanz ++;
	if (myanz > anz)
	{
	    sprintf(buffer, "LinList::check: anz (%d != %d) passt nicht", myanz, anz);
	    error(buffer, false);
	    return;
	}
    }

    if (old_f->next != NULL)
    {
	error("LinList::check: Listenende.next ungleich NULL", false);
	return;
    }

    if (last != old_f)
    {
	error("LinList::check: last ungleich Listenende", false);
	return;
    }

    if (myanz != anz)
    {
	sprintf(buffer, "LinList::check: anz (%d != %d) passt nicht", myanz, anz);
	error(buffer, false);
    }
}


