// Utilities.h: Miscellaneous utility programs for LOsim

int PackTime(int year, int month, int day, int hr, int min, int sec, int msec, long long *time);
void UnpackTime(long long time, int *year, int *month, int *day, int *hr, int *min, int *sec, int *msec);
int AddTime(double incrementSec, int *year, int *month, int *day, int *hr, int *min, int *sec, int *msec);
double DiffTime(long long time0, long long time1);
int DayNumber(int year, int month, int day);
int GetLine(FILE *inFile, char lineIn[], char line[], int maxLine);
void MultipleBlankRemove(char s[]);
void LowerCase(char s[]);
int RemoveWhitespace(char line[]);
int PatternMatch(char string[], char pattern[]);
int FieldMatches(char field[], const char name[]);
int StringValue(char var[], const char value[], int maxLen);
int DateAndTime(long long *pDateTimeVar, const char value[]);
int IntegerValue(long *pVar, const long loLimit, const long hiLimit, const char value[]);
int DoubleValue(double *pVar, const double loLimit, const double hiLimit, const char value[]);
void CheckStr(char field[], char fieldname[]);
void CheckInt(long field, char fieldname[]);
void CheckDbl(double field, char fieldname[]);
void ErrorPrint(char message[]);

