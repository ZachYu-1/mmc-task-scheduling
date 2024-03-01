#include <iostream>
#include <list>
#include <vector>
#include <algorithm>
using namespace std;

#define num_C 3

template <class T>
T sumArr(T data[], int l)
{
    T sum = 0;
    for(int i = 0; i < l; i++)
    {
        sum += data[i];
    }
    return sum;
}

template <class T>
T maxArr(T data[], int l)
{
    T max = 0;
    for(int i = 0; i < l; i++)
    {
        if(data[i] > max) max = data[i];
    }
    return max;
}

template <class T>
T minArr(T data[], int l)
{
    T min = INT_MAX;
    for(int i = 0; i < l; i++)
    {
        if(data[i] < min) min = data[i];
    }
    return min;
}

template <class T>
T avgArr(T data[], int l)
{
    T sum = sumArr(data, l);
    return sum / l;
}

struct task
{
    int id;
    bool isCloud;
    int localCore;
    double exeTime[num_C];                                                   //execution time on local cores
    int rTimeLocal, fTimeLocal, rTimeWS, fTimeWS, fTimeC, fTimeWR;           //ready time, finish time
    int sTimeLocal, sTimeWS;                                                 //start time
    double weight;
    vector<task*> preList;
    vector<task*> sucList;
    int ready1, ready2;
    bool sched;
    task() : id() {}
    task(int i) : id(i), isCloud(false), rTimeLocal(0), fTimeLocal(0), rTimeWS(0), 
    fTimeWS(0), fTimeC(0), fTimeWR(0), sTimeLocal(0),  sTimeWS(0), ready1(0), ready2(0), sched(false) {}
};

struct sched
{
    int id;
    int vTar;
    int kOri;
    int kTar;
    double schedEng;
    int schedTime;
    double engRatio;
    bool tLimit;
    vector<int> schedList[num_C + 1];
    sched(int i, int v, int ko, int kt) : id(i), vTar(v), kOri(ko), kTar(kt), tLimit(false) {}
};

class task_graph
{
    private:
        int numTask;
        int cloudTask = 0;
        int wsTime = 3, cTime = 1, reTime = 1;
        int wTime = wsTime + cTime + reTime;
        int *avaTime;
        double *priList;
        task *tList;
        vector<int>* sList;
    
    public:
        task_graph(int tasks[], int numT, int A[][num_C])
        {
            numTask = numT;
            priList = new double[numT];
            tList = new task[numT];
            sList = new vector<int>[num_C + 1];
            avaTime = new int[num_C + 1];                   //num of cores + 1 ws
            for(int i = 0; i < numT; i++)
            {
                task newT(tasks[i]);
                for(int j = 0; j < num_C; j++)
                {
                    newT.exeTime[j] = A[i][j];
                }
                tList[i] = newT;
            }
        }

        ~task_graph() 
        {
            delete [] tList;
            delete [] priList;
            delete [] avaTime;
            delete [] sList;
        }

        void addEdge(int srcTask, int desTask)
        {
            tList[srcTask].sucList.push_back(&tList[desTask]);
            tList[desTask].preList.push_back(&tList[srcTask]);
        }

        void calReadyTime(int data, bool cloud)
        {
            int temp_rTime = 0;
            int preSize = tList[data].preList.size();

            if(cloud == true){
                if(preSize == 0) tList[data].rTimeWS = 0;
                else
                {
                    for(const task* preTask : tList[data].preList)
                    {
                        if(temp_rTime < preTask->fTimeLocal && preTask->isCloud == false) temp_rTime = preTask->fTimeLocal;
                        if(temp_rTime < preTask->fTimeWS && preTask->isCloud == true) temp_rTime = preTask->fTimeWS;
                    }
                    tList[data].rTimeWS = temp_rTime;
                }
            }
            else{
                if(preSize == 0) 
                {
                    tList[data].rTimeLocal = 0;
                    tList[data].sTimeLocal = 0;
                }
                else
                {
                    for(const task* preTask : tList[data].preList)
                    {
                        if(temp_rTime < preTask->fTimeLocal && preTask->isCloud == false) temp_rTime = preTask->fTimeLocal;
                        if(temp_rTime < preTask->fTimeWR && preTask->isCloud == true) temp_rTime = preTask->fTimeWR;
                    }
                    tList[data].rTimeLocal = temp_rTime;
                }
            }
        }

        int iniSchedule()
        {
            int i, j, k;
            int totalTime = 0;

            //primary assigment
            for(i = 0; i < numTask; i++)
            {
                int minT;
                minT = minArr(tList[i].exeTime, num_C);
                if(minT > wTime)
                {
                    tList[i].isCloud = true;
                    minT = wTime;
                }
            }

            //task prioritizing
            for(i = 0; i < numTask; i++)
            {
                if(tList[i].isCloud == true) tList[i].weight = wTime;
                else tList[i].weight = avgArr(tList[i].exeTime, num_C);
            }

            for(i = numTask - 1; i >= 0; i--)
            {
                int sucSize = tList[i].sucList.size();
                if(sucSize == 0) priList[i] = tList[i].weight;
                else{             
                    double sucPri[sucSize];
                    j = 0;
                    for(const task* sucTask: tList[i].sucList)
                    {
                        sucPri[j] = priList[sucTask->id];
                        j++;
                    }
                    priList[i] = tList[i].weight + maxArr(sucPri, j + 1);
                }
            }

            //exe unit selection
            for(i = 0; i <= num_C; i++)
            {
                avaTime[i] = 0;
            }
            
            for(i = 0; i < numTask; i++)
            {
                double maxPri = -1;
                for(k = 0; k < numTask; k++)
                {
                    if(priList[k] > maxPri) 
                    {
                        maxPri = priList[k];
                        j = k;
                    }
                }
                priList[j] = -1;
                
                //if selected is cloud task
                if(tList[j].isCloud == true)
                {
                    //calculate ready and start time for ws
                    calReadyTime(j, true);
                    tList[j].sTimeWS = max(avaTime[num_C], tList[j].rTimeWS);
                    cloudTask += 1;
                    //update ws, c, wr finish time of task, ws available time and total time
                    tList[j].fTimeWS = tList[j].sTimeWS + wsTime;
                    tList[j].fTimeC = tList[j].fTimeWS + cTime;
                    tList[j].fTimeWR = tList[j].fTimeC + reTime;
                    avaTime[num_C] = tList[j].fTimeWS;
                    totalTime = maxArr(avaTime, num_C + 1);
                    int wsFinishTime = avaTime[num_C] + cloudTask*(cTime + reTime);
                    if(totalTime < wsFinishTime) totalTime = wsFinishTime;
                }
                //if selected is local task
                else
                {
                    calReadyTime(j, true);
                    calReadyTime(j, false);
                    //est time on each core and cloud
                    int temp_fTime[num_C + 1] = {0};
                    for(k = 0; k < num_C; k++)
                    {
                        tList[j].localCore = k;
                        tList[j].sTimeLocal = max(tList[j].rTimeLocal, avaTime[k]);
                        temp_fTime[k] = tList[j].sTimeLocal + tList[j].exeTime[k];
                    }

                    tList[j].sTimeWS = max(avaTime[num_C], tList[j].rTimeWS);
                    temp_fTime[num_C] = tList[j].sTimeWS + wTime;
                    totalTime = minArr(temp_fTime, num_C + 1);
                    k = 0;
                    while(temp_fTime[k] != totalTime) k++;

                    //schedule on local
                    if(k != num_C)
                    {
                        tList[j].localCore = k;
                        tList[j].sTimeLocal = max(tList[j].rTimeLocal, avaTime[k]);
                        //update local finish time of task, local core available time
                        tList[j].fTimeLocal = tList[j].sTimeLocal + tList[j].exeTime[k];
                        avaTime[k] = tList[j].fTimeLocal;
                    }
                    //or reschedule on cloud
                    else if(k == num_C)
                    {
                        tList[j].isCloud = true;
                        tList[j].sTimeWS = max(avaTime[num_C], tList[j].rTimeWS);
                        cloudTask += 1;
                        //update ws, c, wr finish time of task, ws available time
                        tList[j].fTimeWS = tList[j].sTimeWS + wsTime;
                        tList[j].fTimeC = tList[j].fTimeWS + cTime;
                        tList[j].fTimeWR = tList[j].fTimeC + reTime;
                        avaTime[num_C] = tList[j].fTimeWS;
                    }
                    
                }
            }

            return totalTime;
        }

        void recalculateTime(vector<int> s[num_C + 1])
        {
            int i, j, k;
            for(i = 0; i <= num_C; i++) avaTime[i] = 0;
            for(i = 0; i < numTask; i++)
            {
                if(tList[i].isCloud == true) tList[i].weight = wTime;
                else tList[i].weight = avgArr(tList[i].exeTime, num_C);
            }
            for(i = numTask - 1; i >= 0; i--)
            {
                int sucSize = tList[i].sucList.size();
                if(sucSize == 0) priList[i] = tList[i].weight;
                else{             
                    double sucPri[sucSize];
                    j = 0;
                    for(const task* sucTask: tList[i].sucList)
                    {
                        sucPri[j] = priList[sucTask->id];
                        j++;
                    }
                    priList[i] = tList[i].weight + maxArr(sucPri, j + 1);
                }
                tList[i].sTimeLocal = 0;
                tList[i].sTimeWS = 0;
                tList[i].rTimeLocal = 0;
                tList[i].rTimeWS = 0;
                tList[i].fTimeLocal = 0;
                tList[i].fTimeWR = 0;
                tList[i].fTimeWS = 0;
                tList[i].fTimeC = 0;
            }
            for(i = 0; i < numTask; i++)
            {
                double maxPri = -1;
                for(k = 0; k < numTask; k++)
                {
                    if(priList[k] > maxPri) 
                    {
                        maxPri = priList[k];
                        j = k;
                    }
                }
                priList[j] = -1;
                int currK = findPosition(j, s);
                if(currK == num_C)
                {
                    tList[j].isCloud = true;
                    calReadyTime(j, true);
                    tList[j].sTimeWS = max(avaTime[num_C], tList[j].rTimeWS);
                    tList[j].fTimeWS = tList[j].sTimeWS + wsTime;
                    tList[j].fTimeC = tList[j].fTimeWS + cTime;
                    tList[j].fTimeWR = tList[j].fTimeC + reTime;
                    avaTime[num_C] += wsTime;
                }
                else
                {
                    tList[j].isCloud = false;
                    tList[j].localCore = currK;
                    calReadyTime(j, false);
                    tList[j].sTimeLocal = max(tList[j].rTimeLocal, avaTime[currK]);
                    tList[j].fTimeLocal = tList[j].sTimeLocal + tList[j].exeTime[currK];
                    avaTime[currK] = tList[j].fTimeLocal;
                }
            }
        }

        double calculateEng(double p[])
        {
            double totalEng = 0;
            for(int i = 0; i < numTask; i++)
            {
                if(tList[i].isCloud == false) totalEng += p[tList[i].localCore] * (tList[i].exeTime[tList[i].localCore]);
                else totalEng += p[num_C] * wsTime;
            }
            return totalEng;
        }

        void insertTask(task t, vector<int>& s, bool cloud)
        {
            int tIndex = t.id;
            if(cloud == false)
            {
                int rTimetar = t.rTimeLocal;
                int j = 0;
                if(s.size() == 0) s.push_back(tIndex);
                else
                {
                    for(const int& index : s)
                    {
                        if(tList[index].sTimeLocal < rTimetar) j++;
                    }
                    auto it = s.begin();
                    s.insert(it + j, tIndex);
                }
            }
            else
            {
                int rTimetar = t.rTimeWS;
                int j = 0;
                if(s.size() == 0) s.push_back(tIndex);
                else
                {
                    for(const int& index : s)
                    {
                        if(tList[index].sTimeWS < rTimetar) j++;
                    }
                    auto it = s.begin();
                    s.insert(it + j, tIndex);
                }
            }
        }

        void updateReadys(vector<int> s[num_C + 1], bool update, int t)
        {
            if(!update)
            {
                for(int k = 0; k < numTask; k++)
                {
                    for(const task* preTask : tList[k].preList)
                    {
                        if(preTask->sched == false) tList[k].ready1 += 1;
                    }
                }
            }
            else
            {
                for(const task* sucTask : tList[t].sucList)
                {
                    tList[sucTask->id].ready1 -= 1;
                }
            }
            
            for(int k = 0; k <= num_C; k++)
            {
                vector<int> ns = s[k];
                bool flag = true;
                for(const int& index: ns)
                {
                    if(flag == false) tList[index].ready2 = 1;
                    else
                    {
                        tList[index].ready2 = 0;
                        if(tList[index].sched == false) flag = false;
                    }
                }
            }
        }

        int findPosition(int t, vector<int> s[num_C + 1])
        {
            for(int k = 0; k <= num_C; k++)
            {
                vector<int> ns = s[k];
                for(const int& item: ns)
                {
                    if(t == item) return k;
                }
            }
            return -1;
        }

        double taskMigration(int t, double p[])
        {
            int i, j ,k, numMove, numIter;
            int iniTime = t;
            double iniEng = calculateEng(p);
            int maxTime = t * 1.5;
            int totalTime;
            double totalEng = calculateEng(p);
            int tempTime = 0;
            double tempEng = 0;
            bool flagT, flagMax, flagEng;
            vector<int> schStack;
            vector<sched> schList;
            vector<int> newsList[num_C + 1];
            vector<int> tempsList[num_C + 1];

            //construct ini task sequences sList for each core & cloud
            for(i = numTask - 1; i >= 0; i--)
            {
                if(tList[i].isCloud == true) insertTask(tList[i], sList[num_C], true);
                else insertTask(tList[i], sList[tList[i].localCore], false);
            }

            for(i = 0; i < num_C; i++)
            {
                cout << "The ini schedule on Core " << i + 1 <<": " << endl;
                for(const int& index: sList[i])
                {
                    cout << tList[index].id + 1 <<", " << tList[index].sTimeLocal <<", " << tList[index].fTimeLocal << endl;
                }
            }
            cout << "The ini schedule on cloud" << ": " << endl;
            for(const int& index: sList[num_C])
            {
                cout << tList[index].id + 1 <<", " << tList[index].sTimeWS <<", " << tList[index].fTimeWR << endl;
            }
            cout <<"Ini total time: " << iniTime << endl;
            cout <<"Ini total energy : " << iniEng << endl;
            cout << endl;

            numIter = 0;
            flagMax = true;
            while(flagMax)
            {
                //update flags and enter iteration of outer loop
                flagT = false;
                flagMax = false;
                flagEng = false;
                schList.clear();
                numMove = 0;             
                
                //outer loop
                for(i = 0; i < numTask; i++)
                {
                    if(tList[i].isCloud == false)                   //n prime of tasks
                    {
                        for(k = 0; k <= num_C; k++) newsList[k] = sList[k];
                        recalculateTime(sList);
                        int currCore = tList[i].localCore;
                        int oriCore = currCore;
                        for(j = 0; j <= num_C; j++)                 //K cores
                        {
                            if(j != oriCore)
                    //outer loop, Vtar = tList[i], Ktar = j, Kori(1) = currCore
                            {
                                //record this move
                                sched newSched(numMove, i, oriCore, j);
                                
                                //build new sequence after moving vtar
                                k = 0;
                                for(const int& index : newsList[currCore])
                                {
                                    if(tList[i].id != index) k++;
                                    else break;
                                }
                                auto it = newsList[currCore].begin();
                                newsList[currCore].erase(it + k);
                                if(j == num_C) 
                                {
                                    calReadyTime(i, true);
                                    insertTask(tList[i], newsList[j], true);
                                }
                                else 
                                {
                                    calReadyTime(i, false);
                                    insertTask(tList[i], newsList[j], false);
                                }
                                for(k = 0; k <= num_C; k++) newSched.schedList[k] = newsList[k];
                                schList.push_back(newSched);

                                //ini ready1 and ready2, schStack, avaTime
                                for(k = 0; k < numTask; k++) 
                                {
                                    tList[k].sched = false;
                                    tList[k].ready1 = 0;
                                    tList[k].ready2 = 0;
                                }
                                updateReadys(newsList, false, -1);
                                for(k = 0; k < numTask; k++)
                                {
                                    if(tList[k].ready1 == 0 && tList[k].ready2 == 0) schStack.push_back(tList[k].id);
                                }
                                for(k = 0; k <= num_C; k++) avaTime[k] = 0;

                                //schedule all tasks
                                while(schStack.size() > 0)
                                {                              
                                    //pop one task
                                    int index = schStack.back();
                                    schStack.pop_back();
                                    tList[index].sched = true;
                                    int currK = findPosition(index, newsList);
                                    //if it is scheduled on cloud
                                    if(currK == num_C)
                                    {
                                        tList[index].isCloud = true;
                                        calReadyTime(index, true);
                                        tList[index].sTimeWS = max(avaTime[num_C], tList[index].rTimeWS);
                                        tList[index].fTimeWS = tList[index].sTimeWS + wsTime;
                                        tList[index].fTimeC = tList[index].fTimeWS + cTime;
                                        tList[index].fTimeWR = tList[index].fTimeC + reTime;
                                        tList[index].sTimeLocal = 0;
                                        tList[index].fTimeLocal = 0;
                                        avaTime[num_C] = tList[index].fTimeWS;
                                    }
                                    //if it is scheduled on local
                                    else
                                    {
                                        tList[index].isCloud = false;
                                        tList[index].localCore = currK;
                                        calReadyTime(index, false);
                                        tList[index].sTimeLocal = max(tList[index].rTimeLocal, avaTime[currK]);
                                        tList[index].fTimeLocal = tList[index].sTimeLocal + tList[index].exeTime[currK];
                                        tList[index].sTimeWS = 0;
                                        tList[index].fTimeWS = 0;
                                        tList[index].fTimeC = 0;
                                        tList[index].fTimeWR = 0;
                                        avaTime[currK] = tList[index].fTimeLocal;
                                    }
                                    //update ready1 and ready2, schStack
                                    updateReadys(newsList, true, index);
                                    for(k = 0; k < numTask; k++)
                                    {
                                        auto it = find(schStack.begin(), schStack.end(), tList[k].id);
                                        if(it == schStack.end())
                                        {
                                        if(tList[k].ready1 == 0 && tList[k].ready2 == 0 && !tList[k].sched) schStack.push_back(tList[k].id);
                                        }
                                    }
                                }

                                avaTime[num_C] += (cTime + reTime);
                                tempTime = maxArr(avaTime, num_C + 1);
                                tempEng = calculateEng(p);
                                schList[numMove].schedTime = tempTime;
                                schList[numMove].schedEng = tempEng;
                                schList[numMove].engRatio = (iniEng - tempEng)/(tempTime - iniTime);

                                if(tempTime <= iniTime)
                                {
                                    flagT = true;
                                    schList[numMove].tLimit = true;
                                    if(tempEng < totalEng)
                                    {
                                        flagEng = true;
                                        totalEng = tempEng;
                                        totalTime = tempTime;
                                        for(k = 0; k <= num_C; k++) tempsList[k] = newsList[k];                                   
                                    }
                                }
                                else if(tempTime <= maxTime) 
                                {
                                    flagMax = true;
                                    schList[numMove].tLimit = true;
                                }
                                numMove += 1;
                                currCore = j;
                            }
                        }
                    }
                }

                numIter += 1;
                //case1, no increase total time and reduce energy
                if(flagT && flagEng)
                {
                    for(k = 0; k <= num_C; k++) sList[k] = tempsList[k];
                }
                //case2, increase total time but less than max
                else if((!flagT || (flagT && !flagEng)) && flagMax)
                {
                    double tempRatio = 0;
                    int tempID, tempVtar;
                    bool tempFLag = false;
                    for(const sched s: schList)
                    {
                        if((tempRatio < s.engRatio) && s.tLimit && (s.schedEng < totalEng))
                        {
                            tempRatio = s.engRatio;
                            tempID = s.id;
                            tempFLag = true;
                            tempVtar = s.vTar;
                        }
                    }
                    if(tempFLag)
                    {
                        for(k = 0; k <= num_C; k++) sList[k] = schList[tempID].schedList[k];
                        totalTime = schList[tempID].schedTime;
                        totalEng = schList[tempID].schedEng;
                    }
                    else flagMax = false;
                    //case3, total time larger than max, end while
                }

                recalculateTime(sList);
            
                cout << "Number of iteration: " << numIter << endl;
                for(i = 0; i < num_C; i++)
                {
                    cout << "The tasks scheduled on Core " << i + 1 <<": " << endl;
                    for(const int& index: sList[i])
                    {
                        cout << tList[index].id + 1 <<", " << tList[index].sTimeLocal <<", " << tList[index].fTimeLocal << endl;
                    }
                }
                cout << "The tasks scheduled on cloud: " << endl;
                for(const int& index: sList[num_C])
                {
                    cout << tList[index].id + 1 <<", " << tList[index].sTimeWS <<", " << tList[index].fTimeWR << endl;
                }
                cout <<"Total time after iteration " << numIter << " : " << totalTime << endl;
                cout <<"Total energy after iteration " << numIter << " : " << totalEng << endl;
                cout << endl;
            }
            return totalEng;
        }
};

int main()
{
    int iniTotalTime;
    int tasks[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    double pList[] = {1, 2, 4, 0.5};
    int core_matrix[10][3] = {
    {9, 7, 5},
    {8, 6, 5},
    {6, 5, 4},
    {7, 5, 3},
    {5, 4, 2},
    {7, 6, 4},
    {8, 5, 3},
    {6, 4, 2},
    {5, 3, 2},
    {7, 4, 2}};

    task_graph newGraph(tasks, 10, core_matrix);
    newGraph.addEdge(0, 1);
    newGraph.addEdge(0, 2);
    newGraph.addEdge(0, 3);
    newGraph.addEdge(0, 4);
    newGraph.addEdge(0, 5);
    newGraph.addEdge(1, 7);
    newGraph.addEdge(1, 8);
    newGraph.addEdge(2, 6);
    newGraph.addEdge(3, 7);
    newGraph.addEdge(3, 8);
    newGraph.addEdge(4, 8);
    newGraph.addEdge(5, 7);
    newGraph.addEdge(6, 9);
    newGraph.addEdge(7, 9);
    newGraph.addEdge(8, 9);

    iniTotalTime = newGraph.iniSchedule();
    newGraph.taskMigration(iniTotalTime, pList);

}