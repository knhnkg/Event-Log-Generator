/* Program to discover a process model from an event log.

  Skeleton program written by Artem Polyvyanyy, artem.polyvyanyy@unimelb.edu.au,
  August 2022, with the intention that it be modified by students
  to add functionality, as required by the assignment specification.

*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>

/*Global variable definition*************************************/

#define TRACE_END 1
#define REPEATED 2
#define UNIQUE 3
#define CODE 256
#define NONE 4

/*Type definition************************************************/

typedef unsigned int action_t;  // an action is identified by an integer

typedef struct event event_t;   // an event ...
struct event {                  // ... is composed of ...
    action_t actn;              // ... an action that triggered it and ...
    event_t* next;              // ... a pointer to the next event in the trace
};

typedef struct {                // a trace is a linked list of events
    event_t* head;              // a pointer to the first event in this trace
    event_t* foot;              // a pointer to the last event in this trace
    int      freq;              // the number of times this trace was observed
} trace_t;

typedef struct {                // an event log is an array of distinct traces
                                //     sorted lexicographically
    trace_t **trcs;             // an array of traces
    int      ndtr;              // the number of distinct traces in this log
    int      cpct;              // the capacity of this event log as the number
                                //     of  distinct traces it can hold
} log_t;

typedef struct {                // a total trace struct is composed of array 
    trace_t **trcs;             // of pointers to all traces identified from
                                // the file (includes repeated ones), 
    int      nttr;              // the total number of traces and 
    int      cpct;              // capacity = total num of traces it can hold
} total_trace_t;

typedef struct {                // sequence table is composed of an array of 
    action_t **A;               // pointers to array A, that forms a table and
    int       e_num;            // the number of distinct events. 
} table_t;

typedef action_t** DF_t;        // a directly follows relation over actions

/*Function declaration********************************************/

int get_action(); 
trace_t *build_trace(int freq[]);
trace_t *create_empty_trace(void); 
trace_t *insert_at_foot(trace_t *trace, int value); 
void free_trace(trace_t *trace); 
int print_trace (trace_t *t); 
void log_init(log_t *L);
void total_trace_init(total_trace_t *L);
int distinct_trace_check(trace_t *t, log_t log, int i);
int trace_compare (trace_t *t, trace_t *a); 
void most_freq_trace(log_t log, int i);
void event_frequency(int freq[]);
int distinct_event(int freq[]);
int total_event(int freq[]);
int qcmp(const void *aa, const void *bb);

int w_calc(int pd, int sup1, int sup2);
int pd_calc(int sup1, int sup2);
int max_calc(int sup1, int sup2);

void table_init(table_t *M, int distinct_event_num);
void table_col_row(table_t table, int freq[]);
int table_event(int freq[], int k);
int table_support (trace_t *t, table_t table, int code);
void support_freq (table_t table, int actn1, int actn2, int code);
void event_freq_removed(int freq[], int event1, int event2);
void total_trcs_replace (trace_t *t, int freq[], int sup1, int sup2, int code);
int candidate_checker(table_t table);

int support_build (total_trace_t total_trcs, int freq[], int abstract[], 
                   int distinct_event_num, int total_trace_counter, 
                   int code, int Q);
    
/****************************************************************/

/* Main program controls all actions of the program. */
int main(int argc, char *argv[]) {

    // STAGE 0:
    printf("==STAGE 0============================\n");
  
    trace_t *t;
    int total_trace_counter = 0;
    int i = 0;
    
    // This struct contains infor for all DISTINCT traces. 
    log_t log;
    log_init(&log);
    
    // This struct contains infor for all traces (includes repeated ones).  
    total_trace_t total_trcs;
    total_trace_init(&total_trcs);
    
    // This array records the frequency of each event. 
    int freq[512] = {0};
    
    while ( ((t = build_trace(freq)) != NULL) ) {
        // Another trace identified, increase the counter by 1.
        total_trace_counter += 1;
        
        // Check if the struct capacity is full and increase if needed.
        if (total_trcs.nttr == total_trcs.cpct) {
            total_trcs.cpct *= 2;
            total_trcs.trcs = (trace_t**)realloc(total_trcs.trcs, 
                               total_trcs.cpct * sizeof(*(total_trcs.trcs)));
            assert(total_trcs.trcs);
        }
        // Insert identifed trace into struct. 
        total_trcs.trcs[total_trcs.nttr] = t;
        total_trcs.nttr++;
        
        // When log struct is empty, it indicates that this identified struct
        // is definetely DISINTCT, so add to the struct without comparison.
        if (log.ndtr == 0) {
            log.trcs[log.ndtr] = t;
            log.trcs[log.ndtr]->freq++;
            log.ndtr++;
        }
        
        // Compare the identified trace with all other distinct traces in the
        // log struct, before adding it.
        else {
            int check = distinct_trace_check(t, log, i);
            
            // If repeated, do nothing about the trace. 
            if (check == REPEATED) {
            }
            
            // If DISTINCT, check if the capacity of the struct is full and if
            // needed, increase the capacity first. Then add the trace. 
            else if (check == UNIQUE) {
                if (log.ndtr==log.cpct) {
                    log.cpct *= 2;
                    log.trcs = (trace_t**)realloc(log.trcs, 
                                log.cpct * sizeof(*(log.trcs)));
                    assert(log.trcs);
                } 
                log.trcs[log.ndtr]= t; 
                log.trcs[log.ndtr]->freq++;
                log.ndtr++;
                i += 1;
            }
        }      
    } 
    
    qsort(log.trcs, log.ndtr, sizeof( * (log.trcs)), qcmp); // Sort log.    
    int distinct_event_num = distinct_event(freq);
    printf("Number of distinct events: %d\n", distinct_event_num);    
    printf("Number of distinct traces: %d\n", i+1);   
    int total_event_num = total_event(freq);
    printf("Total number of events: %d\n", total_event_num);    
    printf("Total number of traces: %d\n", total_trace_counter);   
    most_freq_trace(log,i);   
    event_frequency(freq);
    
    // STAGE 1: 
    printf("==STAGE 1============================\n");
    
    
    int Q = 0;
    int abstract[300] = {0}; // Keeps record of abstract sequence values. 
    while (Q <= 1) {
        support_build (total_trcs, freq, abstract, distinct_event_num, 
                       total_trace_counter, CODE, Q);  
        printf("=====================================\n");
        Q += 1;
    }
    
    return 0;
}

/* Initialise the members of "log_t" typedef struct. */
void log_init(log_t *L) {
    L->cpct = 10; // Can be any small beginning value for capacity.
    L->ndtr = 0;  // Begin from 0 to compare with cpct.
    L->trcs = malloc(L->cpct * sizeof(*(L->trcs)));
    assert(L->trcs);
}

/* Initialise the members of "total_trace_t" typedef struct. */
void total_trace_init(total_trace_t *L) {
    L->cpct = 10; // Can be any small beginning value for capacity.
    L->nttr = 0;  // Begin from 0 to compare with cpct.
    L->trcs = malloc(L->cpct * sizeof(*(L->trcs)));
    assert(L->trcs);
}

/* Compares if two argument traces are same or different. */
int trace_compare(trace_t *t, trace_t *a) {
    event_t *p, *q;
    p = t->head; // Begin from head for both traces. 
    q = a->head;    
    while (p != NULL && q != NULL) {
        // Difference detected, stop the loop. 
        if (p->actn < q->actn) return -1;
        // Difference detected, stop the loop. 
        else if (p->actn > q->actn) return 1;
        // Identical, move to next event for both traces. 
        p = p->next;
        q = q->next;
    }
    // All events of two traces are identical. 
    if (p == NULL && q == NULL) return 0;
    // One trace is longer than the other, hence two are not identical.
    if (p == NULL) return -1;
    // One trace is longer than the other, hence two are not identical.
    return 1;
}

/* Compares input traces. */
int qcmp(const void *aa, const void *bb) {
    trace_t *a= (trace_t *) *((trace_t **) aa);
    trace_t *b= (trace_t *) *((trace_t **)bb);  
    if (a->freq>b->freq) return -1;
    if (a->freq<b->freq) return +1;
    return trace_compare(a, b);
}

/* Determines if the input trace is distinct from all other traces. */
int distinct_trace_check(trace_t *t, log_t log, int i) {
    int j = 0;    
    while (j <= i) {
        trace_t *a = log.trcs[j];
        // Compare input trace to all other traces in the log struct.
        int check = trace_compare(t, a);
        if (check == 0) {
            // If identical to one of them, increase the freq of that trace
            // and do not add the input trace to the log again. 
            log.trcs[j]->freq++;
            return REPEATED;
        }
        j += 1;       
    }
    return UNIQUE;
}

/* Identifies the trace(s) that occured the most and their frequency. */
void most_freq_trace(log_t log, int i) {
    int freq = 0;   
    for (int j=0; j <= i; j++) {
        // Find the max frequency by comparing current freq to max freq. 
        if (log.trcs[j]->freq > freq) {
            // Replace max freq by current freq if its greater. 
            freq = log.trcs[j]->freq;   
        }
    }    
    printf("Most frequent trace frequency: %d\n", freq);    
    for (int j=0; j <= i; j++) {
        // Identify traces that has the max freq identified above.
        if (log.trcs[j]->freq == freq) {
            print_trace(log.trcs[j]);
        }
    }    
}

/* Prints the frequency of occurence of each distinct event. */
void event_frequency(int freq[]) {
    for (int t = 0; t <= 512; t++) {
        // The fact that freq of a event is greater than 0 indicates that
        // this event occured at least once in at least one trace. 
        if ((freq[t] > 0) && (freq[t] < 512)) {
            // print the event char and its freq. 
            printf("%c = %d\n", t, freq[t]);
        }
    }
}

/* Counts the number of distinct events. */
int distinct_event(int freq[]) {
    int distinct_event_counter = 0;    
    for (int t = 0; t <= 512; t++) {
        // The fact that freq of a event is greater than 0 indicates that
        // this event occured at least once in at least one trace. 
        if ((freq[t] > 0) && (freq[t] < 512)) {
            // Increase counter by 1. 
            distinct_event_counter += 1;
        }
    }
    return distinct_event_counter;
}

/* Counts the total number of events. */
int total_event(int freq[]) {
    int total_event_counter = 0;    
    for (int t = 0; t <= 512; t++) { 
        // The fact that freq of a event is greater than 0 indicates that
        // this event occured at least once in at least one trace. 
        if ((freq[t] > 0) && (freq[t] < 512)) {
            // Add freq of this trace to the total counter. 
            total_event_counter = total_event_counter + freq[t];
        }
    }
    return total_event_counter;
}

/* Print trace. */
int print_trace (trace_t *t) {
    event_t *p;
    p = t->head;                  // Begine at the head of trace.
    while(p != NULL) {
        printf("%c", p->actn);    // Move to action and print.
        p = p->next;              // Move to the next action, if not NULL. 
    }
    printf("\n");
    return 0;
}

/* Gather each event extracted from file into corresponding trace struct. */
trace_t *build_trace(int freq[]) {
    trace_t *t= create_empty_trace();           // Create new trace struct.
    int c;
    int num = 0;    
    while ((c=get_action()) && isalpha(c) ) {   // Extract an event. 
        t = insert_at_foot(t, c);               // Add to corresponding trace.
        freq[c] += 1;                           // Increase freq of event. 
        num++;  
    }
    if (num>0) {                                // Trace is filled & needed.
        return t;
    } 
    else {                                      // Trace is empty and so not...
        free_trace(t);                          // ...required, so free it. 
        return NULL;                              
    }
}

/* Create a new trace. */
trace_t *create_empty_trace(void) {
    trace_t *trace;
    trace = (trace_t*)malloc(sizeof(*trace));
	assert(trace!=NULL);
	trace->head = trace->foot = NULL;
	return trace;
}

/* Insert extracted event to the foot(end) of corresponding trace. */
trace_t *insert_at_foot(trace_t *trace, int value) {
	event_t* new;
	new = (event_t*)malloc(sizeof(*new));
	assert(trace!=NULL && new!=NULL);
	new->actn = value;
	new->next = NULL;
	if (trace->foot==NULL) {
		trace->head = trace->foot = new;
	} else {
		trace->foot->next = new;
		trace->foot = new;
	}
	return trace;
}

/* Free trace elements if not needed anymore. */
void free_trace(trace_t *trace) {
	event_t *curr, *prev;
	assert(trace!=NULL);
	curr = trace->head;
	while (curr) {
		prev = curr;
		curr = curr->next;
		free(prev);
	}
	free(trace);
}

/* Extract single event character at a time. */
int get_action() {
    int c;    
    while ((c=getchar()) != EOF && c!='\n' && !isalpha(c))
    if (c==EOF) return EOF;            // Found end of file.
    if (c=='\n') return TRACE_END;     // Found end of a trace.
    return c;                          // Found event character.
}

/* Calculate the weight of sup(u,v) and sup(v,u). */
int w_calc(int pd, int sup1, int sup2) {
    int max = max_calc(sup1, sup2);     // Find the max value between two. 
    int w = abs(50 - pd) * max;         // Formula to calculate w.
    return w;
}

/* Calculate the percent difference of sup(u,v) and sup(v,u). */
int pd_calc(int sup1, int sup2) {
    int max = max_calc(sup1, sup2);         // Find the max value between two. 
    int pd = (100 * abs(sup1 - sup2)) / max;        // Formula to calculate w.
    return pd;
}

/* Identifies the max value out of two input integers. */
int max_calc(int sup1, int sup2) {
    if (sup1 > sup2) return sup1; // LHS is larger than RHS, return LHS. 
    else if (sup1 < sup2) return sup2;  // RHS is larger than LHS, return RHS.
    return 0; // LHS and RHS are identical. 
}

/* Initialise the members of "table_t" typedef struct. */
void table_init(table_t *M, int distinct_event_num) {
    M -> e_num = distinct_event_num;
    M -> A = malloc((M->e_num+1) * sizeof(*(M->A)));
    assert(M->A);
    for (int i = 0; i <= (M->e_num); i++) {
        M->A[i] = malloc((M->e_num+1) * sizeof(*(M->A[i])));
        assert(M->A[i]);
    }
}

/* Inserts distinct event characters into column headers and row headers. */
void table_col_row(table_t table, int freq[]) {
    int k = 0;
    int y = 0;  
    // ROWS (VERTICAL)
    for (int i=0; i <= 0; i++) {
        for (int j=1; j <= table.e_num ; j++) { 
            // Found a distinct event char that occurs at least once. 
            y = table_event(freq, k);
            // Set it to one of the row headers. 
            table.A[i][j] = y;
            k = y + 1;
        }
    }
    
    int l = 0;
    int x = 0;
    // COLUMNS (HORIZONTAL)
    for (int i=1; i <= table.e_num; i++) {
        for (int j=0; j <= 0; j++) {
            // Found a distinct event char that occurs at least once. 
            x = table_event(freq, l);
            // Set it to one of the column headers. 
            table.A[i][j] = x;
            l = x + 1;
        }
    }
}

/* Identifies the distinct events that occured in at least one trace, 
   one at a time. */
int table_event(int freq[], int k) {
    while (k <= 512) {
        if ((freq[k] > 0) && (freq[k] < 512)) {
            return k;
        }
        k ++;
    }
    return 0;
}

/* Finds all sup(u,v) pattern in each trace. */
int table_support (trace_t *t, table_t table, int code) {
    event_t *p;
    int actn1;
    int actn2;
    p = t->head;                 // Begin from the head of the trace.  
    while(p != NULL) {
        actn1 = p->actn;         // Found the first action of a sup(u,v). 
        p = p->next;
        if (p != NULL) {
            actn2 = p->actn;     // Found the following action of a sup(u,v). 
            support_freq(table, actn1, actn2, code);
        }
    }
    return 0;
}

/* Record frequency of each sup(u,v) pattern. */
void support_freq (table_t table, int actn1, int actn2, int code) {
    int row = 0;
    int col = 0;
    
    // ROWS (VERTICAL):
    for (int i=0; i <= 0; i++) {
        for (int j=1; j <= table.e_num ; j++) {
            // Find the table index that corresponds to action input. 
            if ((int)table.A[i][j] == actn1) {
                row = j; // Found the index for the row. 
                break;
            }
        }
    }

    // COLUMNS (HORIZONTAL):
    for (int i=1; i <= table.e_num; i++) {
        for (int j=0; j <= 0; j++) {
            // Find the table index that corresponds to action input.
            if ((int)table.A[i][j] == actn2) {
                col = i; // Found the index for the column.
                break;
            }
        }
    }
    
    // Ignore, if the position is for abstract sequence. 
    if ( (actn1 == code) && (actn2 == code) ) {
        table.A[row][col] = 0;
    }
    
    // Increase the freq at determined table index by 1. 
    else {
        table.A[row][col] += 1;
    }
}

/* Print frequency of all distinct events, except the two events of 
   "most likely sequential pattern". */
void event_freq_removed(int freq[], int event1, int event2) {
    for (int t = 0; t <= 122; t++) {
        // Print freq of all event that occured in at least once in a trace,
        // except for events that equals to one of the action from sup(u,v). 
        if ((freq[t] > 0) && (freq[t] < 512) && 
            (t != event1) && (t != event2)) {
            printf("%c = %d\n", t, freq[t]);
        }
    }
}

/* Replace the two events of "most likely sequential pattern" with their 
   corresponding code (eg.256) and modify their frequencies. */
void
total_trcs_replace(trace_t *t, int freq[], int sup1, int sup2, int code) {
    event_t *p;
    int actn;
    p = t->head; // Begin from the head of the trace. 
   
    while(p != NULL) {
        actn = p->actn;
        if ( (actn == sup1) || (actn == sup2) ) {
            // If current action equals to one of the action from sup(u,v), 
            // move it's freq to corresponding code and replace it's freq 
            // with -1 for future identification. 
            freq[code] = freq[code] + freq[actn];
            p->actn = code;
            freq[actn] = -1;
        }
        p = p->next;  // Move to next to check for next action. 
    }
}

/* Check if there is any candidate for "most likely sequential pattern". */
int candidate_checker(table_t table) {
    int sup1=0, sup2=0;
    int max_w = 0;    
    for (int r=1; r <= table.e_num; r++) {
        for (int c=1; c <= table.e_num; c++) {
            // Check for CONDITION 1. For sup(u,v), u != v.
            if (c != r) {
                if ( ((int)table.A[r][0] < 256) && ((int)table.A[0][c]) ) {
                    // Check for CONDITION 2. sup(u,v) > sup(v,u). 
                    if ((int)table.A[r][c] > (int)table.A[c][r]) {
                        int pd = pd_calc((int)table.A[r][c], 
                                         (int)table.A[c][r]);
                        int w = w_calc(pd, (int)table.A[r][c], 
                                           (int)table.A[c][r]);
                        // Check for CONDITION 3. pd(u,v) > 70. 
                        // And to only select one, select one with max weight. 
                        if ((pd > 70) && (w > max_w)) {
                            sup1 = (int)table.A[r][0];
                            sup2 = (int)table.A[0][c];
                            max_w = w;
                        }
                    }
                }
            }
        }
    }
    // If sup1, sup2 and max_w has not been changed from their initial value,
    // this indicates that there is no more candidate. 
    if ( (sup1 == 0) && (sup2 == 0) && (max_w == 0) ) {
        return 0;
    }
    // Otherwise, return 1 to indicate that another candidate exists. 
    return 1; 
}

int support_build (total_trace_t total_trcs, int freq[], int abstract[], 
               int distinct_event_num, int total_trace_counter, 
               int code, int Q) {
    // Initialise support table.
    table_t table; 
    table_init(&table, distinct_event_num-Q);  
    
    // Set all row and column values to 0 first.
    for (int t=0; t <= table.e_num; t++) {
        for (int s=0; s <= table.e_num; s++) {
            table.A[t][s] = 0;
        }
    }    
    // Insert column and row headers with event characters.
    table_col_row(table, freq); 
    
    // Record frequency of each support sequence. 
    for (int z=0; z<total_trace_counter; z++) {
        table_support(total_trcs.trcs[z], table, Q);
    }     
    // Assign blank space to A[0][0].
    table.A[0][0] = 32; 
    
    // Check if there is a candidate that satisfies all conditions. 
    int candidate_check = candidate_checker(table);
    
    if (candidate_check == 1) {
        // Print the table. 
        for (int t=0; t <= table.e_num; t++) {
            for (int s=0; s <= table.e_num; s++) {
                if (((int)table.A[t][s] > total_trace_counter) && 
                    ((int)table.A[t][s] != 255+Q)) {
                    printf("  %3c", table.A[t][s]);
                }
                else {
                    printf("  %3d", table.A[t][s]);
                }
            }
            printf("\n");
        }
        printf("-------------------------------------\n"); 

        // Find the corresponding "most likely sequential pattern". 
        int sup1, sup2, sup_freq, max_w = 0;    
        for (int r=1; r <= table.e_num; r++) {
            for (int c=1; c <= table.e_num; c++) {
                // Check for CONDITION 1. For sup(u,v), u != v.
                if (c != r) {
                    // Check for CONDITION 2. sup(u,v) > sup(v,u). 
                    if ((int)table.A[r][c] > (int)table.A[c][r]) {
                        int pd = pd_calc((int)table.A[r][c], 
                                         (int)table.A[c][r]);
                        int w = w_calc(pd, (int)table.A[r][c], 
                                           (int)table.A[c][r]);
                        // Check for CONDITION 3. pd(u,v) > 70. 
                        // And to only select one, select one with max weight. 
                        if ((pd > 70) && (w > max_w)) {
                            sup1 = (int)table.A[r][0];
                            sup2 = (int)table.A[0][c];
                            sup_freq = (int)table.A[r][c];
                            max_w = w;
                        }
                    }
                }
            }
        }        
        printf("%d = SEQ(%c,%c)\n", code+Q, sup1, sup2);
        printf("Number of events removed: %d\n", sup_freq);
        event_freq_removed(freq, sup1, sup2);
        
        abstract[code+Q] = sup_freq;
        for (int q = 0; q < 300; q++) {
            if (abstract[q] > 0) {
                printf("%d = %d\n", q, abstract[q]);
            }
        }
        
        // Remove/ exclude the actions of "most likely sequential pattern".
        for (int z=0; z<total_trace_counter; z++) {
            total_trcs_replace(total_trcs.trcs[z], freq, sup1, sup2, code+Q);
        }
    }
    
    else if (candidate_check == 0) {
        // Free all arrays of the table A. 
        for (int i = 0; i < table.e_num; i++) {
          free(table.A[i]);
        }
        free(table.A);
        return NONE;
    }
    
    // Free all arrays of the table A. 
    for (int i = 0; i <= table.e_num; i++) {
      free(table.A[i]);
    }    
    free(table.A);   
    return 0;
}

/* algorithms are fun */

/* ============================================================================
   Functions "create_empty_trace", "insert_at_foot" and "free_trace" were
   modified from the Program written by Alistair Moffat, as an example for 
   the book "Programming, Problem Solving, and Abstraction with C", Pearson
   Custom Books, Sydney, Australia, 2002; revised edition 2012. 
   
   See link http://people.eng.unimelb.edu.au/ammoffat/ppsaa/ for further
   information.
   ========================================================================= */

/* ============================================================================
   Functions print_trace, trace_compare and qcmp were modified from the Program
   covered during tutorial material, by the tutor.  
   ========================================================================= */

/* THE END -------------------------------------------------------------------*/