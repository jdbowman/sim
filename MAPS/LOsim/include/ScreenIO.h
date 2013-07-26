#define MAX_LINES 10000
#define LINE_LENGTH 200

void ScreenInit(void);
void ProcessScrollMessage(int scrollMsg, int param);
void OutputToScreen(char x[], int sameLine);
void SkipLine(void);
