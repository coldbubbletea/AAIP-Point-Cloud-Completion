typedef map<string,string> ConfigType;

void TrimSpace(char* str);
void AddConfigFromFile(ConfigType &cr,const char* filename);
void AddConfigFromCmdLine(ConfigType &cr,int argc,char** argv);
void ListConfig(ConfigType &cr);

float getConfigFloat(const char* key,ConfigType &cr,bool required=true,float _default=0);
int getConfigInt(const char* key,ConfigType &cr,bool required=true,int _default=0);
const char* getConfigStr(const char* key,ConfigType &cr,bool required=true,const char* _default=NULL);

void TrimSpace(char* str) {
	if (str==NULL) return;

	char space[]={'\t','\n','\f','\r',' '};
	int pos=0;
	for (int i=0;i<strlen(str);i++) {
		bool found=false;
		for (int j=0;j<5;j++)
			if (str[i]==space[j]) found=true;

		if (!found) {
			str[pos]=str[i];
			pos++;
		}
	}
	str[pos]='\0';
}

void AddConfigFromFile(ConfigType &cr,const char* filename) {
	const int LINE_LEN=1024;
	char line[LINE_LEN],key[LINE_LEN],value[LINE_LEN];
    char filename2[LINE_LEN];

	ifstream br(filename);
  	if (! br.is_open()) { 
        //printf("Error opening file \"%s\"",filename); 
        sprintf(filename2, "../%s", filename);
        br.clear();
        br.open(filename2);
        if (! br.is_open()) {
            printf("Error opening file at location \"%s\" nor \"%s\"",filename, filename2); 
            exit (1); 
        }
    }

	while (br.getline(line,LINE_LEN)){
		if (strstr(line,"//")!=NULL) continue; // remove comments
		char* chPos=strchr(line,'=');
		if (chPos!=NULL) {
			int pos=((int)(chPos-line))/sizeof(char);
			int keyLen=pos;
			int valueLen=strlen(line)-1-keyLen;
			memcpy(key,&line[0],keyLen);	key[keyLen]='\0';
			memcpy(value,&line[pos+1],valueLen);	value[valueLen]='\0';
			TrimSpace(key);	TrimSpace(value);
			cr[key]=value;
		}
	}
	br.close();
}

void AddConfigFromCmdLine(ConfigType &cr,int argc,char** argv) {
	int i=0;
	while (i<argc) {
		while ((i<argc)&&(argv[i][0]!='-')) i++;	// shortcut condition
		if (i+1<argc) {
			char* key=&(argv[i][1]);
			char* value=argv[i+1];
			TrimSpace(key);	TrimSpace(value);
			cr[key]=value;
			i+=2;
		} else
			return;
	}
}

void ListConfig(ConfigType &cr) {
	ConfigType::iterator p=cr.begin();
	while (p!=cr.end()) {
		printf("%s=%s\n",p->first.c_str(),p->second.c_str());
		p++;
	}
}

float getConfigFloat(const char* key,ConfigType &cr,bool required,float _default) {
	float value=_default;
	if (cr.count(key))
		value=atof(cr[key].c_str());
	else {
		if (required) {
			printf("Config key \"%s\" not found\n",key);
			exit(1);
		}
	}
	return value;
}

int getConfigInt(const char* key,ConfigType &cr,bool required,int _default) {
	int value=_default;
	if (cr.count(key))
		value=atoi(cr[key].c_str());
	else {
		if (required) {
			printf("Config key \"%s\" not found\n",key);
			exit(1);
		}
	}
	return value;
}

const char* getConfigStr(const char* key,ConfigType &cr,bool required,const char* _default) {
	const char* value=_default;
	if (cr.count(key))
		value=cr[key].c_str();
	else {
		if (required) {
			printf("Config key \"%s\" not found\n",key);
			exit(1);
		}
	}
	return value;
}

