/* Continuous Time Bayesian Network Reasoning and Learning Engine
 * Copyright (C) 2009 The Regents of the University of California
 *
 * see docs/AUTHORS for contributor list
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "streamserial.h"

#include "nullptr03.h"
#include <cstring>

namespace ctbn {

// the (few) definitions for that in streamserial.h

std::map<std::string, StreamObj::CreatorFn> &StreamObj::loadtable() {
	static std::map<std::string, StreamObj::CreatorFn> *_loadtable
		= new std::map<std::string, StreamObj::CreatorFn>();
	return *_loadtable;
}

StreamObj::StreamObj() {
}

StreamObj::~StreamObj() {
}

const char *StreamObj::myID() const {
	return "StreamObj";
}

const char *StreamObj::IDname() {
	static const char *ret = StreamObj::addID("StreamObj",StreamObj::Create);
	return ret;
}

void StreamObj::SaveOld(std::ostream &os) const {
}

void StreamObj::SaveOldPtr(std::ostream &os) const {
	os << myID() << os.fill();
	SaveOld(os);
}

StreamObj *StreamObj::LoadOldPtr(std::istream &is) {
	std::string name;
	is >> name;
	std::map<std::string, CreatorFn>::iterator rec
			= loadtable().find(name);
	if (rec == loadtable().end()) return nullptr03;
	else return rec->second(is);
        
}

StreamObj *StreamObj::Create(std::istream &is) {
	return new StreamObj();
}

const char *StreamObj::addID(const char *classname, CreatorFn fn) {
	loadtable()[std::string(classname)] = fn;
	return classname;
}

const char *StreamObj::addID(const char *classname, const char *classname2,
					CreatorFn fn) {
	char *newname = new char[strlen(classname)+strlen(classname2)+3];
	strcpy(newname,classname);
	strcat(newname,"<");
	strcat(newname,classname2);
	strcat(newname,">");
	loadtable()[std::string(newname)] = fn;
	return newname;
}

char *striptypename(const char *name) {
     if (name && name[0]=='N') {
          int i = atoi(name+1);
          int j = 1;
          for(;name[j]<='9' && name[j]>='0';j++)
               ;
		char *ret = new char[strlen(name)-j-i+1];
		i +=j;
		int b = 1,k=0;
		for(;name[i];i++) {
			if (name[i]=='I' && name[i+1]=='N' && name[i+2]=='S'
				&& name[i+3]=='_') {
				ret[k++] = 'I';
				b++;
				i += 3;
			} else ret[k++] = name[i];
		}
		ret[k-b] = 0;
          return ret;
     } else return strcpy(new char[strlen(name)+1],name);
}

const char *StreamObj::_forcename = StreamObj::IDname();

/*
char *StreamObj::killSpaces(const char *classname) {
	char *ret = new char[strlen(classname)+1];
	strcpy(ret,classname);
	for(int i=0;ret[i];i++)
		if (ret[i]==' ') ret[i]='_';
	return ret;
}
*/

} // end of ctbn namespace
