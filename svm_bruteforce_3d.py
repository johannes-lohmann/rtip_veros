import matplotlib.pyplot as pl
import numpy as np
import xarray as xr
import matplotlib as mpl
from sklearn.metrics import accuracy_score
from sklearn import svm, datasets


def MAIN():
        variables = ['salt_sub_NA', 'salt_sub_SA', 'temp_sub_NA', 'temp_sub_SA', 'sst_NA', 'sst_SA', 'sss_NA', 'sss_SA', 'rho_sub_NA', 'rho_sub_SA', 'rho_NA', 'rho_SA', 'seaice', 'AMOC', 'AABW', 'eq', 'delta_rho', 'delta_rho_sub', 'delta_sst', 'delta_temp_sub', 'delta_salt_sub']

        M=60
        amoc=[]; eq=[]; aabw=[]
        edge_reals = [6, 11, 14, 15, 17, 19, 20, 23, 24, 29, 35, 40, 46, 54, 57]
        for i in range(M):
                ot = xr.open_dataarray("overturning_rtipz%s.nc"%i)
                t_amoc0, amoc0, eq0, aabw0 = amoc_aabw_timeseries_ext(ot)
                amoc.append(amoc0);eq.append(eq0);aabw.append(aabw0)
                

        scores = []; count = 0
        best_lin1 = 0.; best_lin2 = 0.; best_rbf = 0.; best_poly = 0.
        for j in range(len(variables)):
                X=np.empty((M,3)); y=np.empty(M)
                print '----------------------------------------------'
                print '----------------------------------------------'
                print '----------------------------------------------'

                for k in range(len(variables)):
                        print '----------------------------------------------'
                        print '----------------------------------------------'
                        for l in range(len(variables)):
                                
                                if (k>j and l!=k and l!=j):
                                        count+=1
                                        print j, k, l
                                        print '----------------------------------------------'
                                        for i in range(M):
                                                st = xr.open_dataset("salt_temp_rtipz%s.nc"%i)
                                                amoc0=amoc[i]
                                                if variables[j]=='eq':
                                                        X0 = eq[i]
                                                elif variables[j]=='AMOC':
                                                        X0 = amoc[i]
                                                elif variables[j]=='AABW':
                                                        X0 = aabw[i]
                                                elif variables[j]=='delta_rho':
                                                        X0 = st['rho_NA'].values-st['rho_SA'].values
                                                elif variables[j]=='delta_sst':
                                                        X0 = st['sst_NA'].values-st['sst_SA'].values
                                                elif variables[j]=='delta_rho_sub':
                                                        X0 = st['rho_sub_NA'].values-st['rho_sub_SA'].values
                                                elif variables[j]=='delta_temp_sub':
                                                        X0 = st['temp_sub_NA'].values-st['temp_sub_SA'].values
                                                elif variables[j]=='delta_salt_sub':
                                                        X0 = st['salt_sub_NA'].values-st['salt_sub_SA'].values
                                                else:
                                                        X0 = st[variables[j]].values

                                                if variables[k]=='eq':
                                                        Y0 = eq[i]
                                                elif variables[k]=='AMOC':
                                                        Y0 = amoc[i]
                                                elif variables[k]=='AABW':
                                                        Y0 = aabw[i]
                                                elif variables[k]=='delta_rho':
                                                        Y0 = st['rho_NA'].values-st['rho_SA'].values
                                                elif variables[k]=='delta_sst':
                                                        Y0 = st['sst_NA'].values-st['sst_SA'].values
                                                elif variables[k]=='delta_rho_sub':
                                                        Y0 = st['rho_sub_NA'].values-st['rho_sub_SA'].values
                                                elif variables[k]=='delta_temp_sub':
                                                        Y0 = st['temp_sub_NA'].values-st['temp_sub_SA'].values
                                                elif variables[k]=='delta_salt_sub':
                                                        Y0 = st['salt_sub_NA'].values-st['salt_sub_SA'].values
                                                else:
                                                        Y0 = st[variables[k]].values

                                                if variables[l]=='eq':
                                                        Z0 = eq[i]
                                                elif variables[l]=='AMOC':
                                                        Z0 = amoc[i]
                                                elif variables[l]=='AABW':
                                                        Z0 = aabw[i]
                                                elif variables[l]=='delta_rho':
                                                        Z0 = st['rho_NA'].values-st['rho_SA'].values
                                                elif variables[l]=='delta_sst':
                                                        Z0 = st['sst_NA'].values-st['sst_SA'].values
                                                elif variables[l]=='delta_rho_sub':
                                                        Z0 = st['rho_sub_NA'].values-st['rho_sub_SA'].values
                                                elif variables[l]=='delta_temp_sub':
                                                        Z0 = st['temp_sub_NA'].values-st['temp_sub_SA'].values
                                                elif variables[l]=='delta_salt_sub':
                                                        Z0 = st['salt_sub_NA'].values-st['salt_sub_SA'].values
                                                else:
                                                        Z0 = st[variables[l]].values

                                                X[i,:] = [X0[14], Y0[14], Z0[14]]

                                                if i in edge_reals:
                                                        y[i] = 1
                                                elif amoc0[-1]<4.:
                                                        y[i] = 2
                                                else:
                                                        y[i] = 0

                                        X[:,0] = (X[:,0]-np.mean(X[:,0]))/np.std(X[:,0])
                                        X[:,1] = (X[:,1]-np.mean(X[:,1]))/np.std(X[:,1])
                                        X[:,2] = (X[:,2]-np.mean(X[:,2]))/np.std(X[:,2])

                                        C = 1.  # SVM regularization parameter
                                        models = (svm.SVC(kernel='linear', C=C),
                                                  svm.LinearSVC(C=C, max_iter=10000),
                                                  svm.SVC(kernel='rbf', gamma=0.7, C=C),
                                                  svm.SVC(kernel='poly', degree=4, gamma='auto', C=C))
                                        models = (clf.fit(X, y) for clf in models)

                                        titles = ('SVC with linear kernel',
                                          'LinearSVC (linear kernel)',
                                          'SVC with RBF kernel',
                                          'SVC with polynomial (degree 4) kernel')
                                        for clf, title in zip(models, titles):
                                                y_pred = clf.predict(X)
                                                accuracy = accuracy_score(y, y_pred)
                                                scores.append(accuracy)
                                                if accuracy>0.65:
                                                        print variables[k], ' - ', variables[j], ' - ', variables[l], ': ', title, ': ', round(accuracy,3)

                                                if title=='SVC with linear kernel':
                                                        if accuracy>best_lin1:
                                                                combo_lin1 = [variables[k], variables[j], variables[l]]
                                                                best_lin1=accuracy

                                                elif title=='LinearSVC (linear kernel)':
                                                        if accuracy>best_lin2:
                                                                combo_lin2 = [variables[k], variables[j], variables[l]]
                                                                best_lin2=accuracy

                                                elif title=='SVC with RBF kernel':
                                                        if accuracy>best_rbf:
                                                                combo_rbf = [variables[k], variables[j], variables[l]]
                                                                best_rbf=accuracy

                                                else:
                                                        if accuracy>best_poly:
                                                                combo_poly = [variables[k], variables[j], variables[l]]
                                                                best_poly=accuracy

        print max(scores)
        print count
        
        print 'Best: SVC with linear kernel'
        print best_lin1
        print combo_lin1

        print 'Best: LinearSVC (linear kernel)'
        print best_lin2
        print combo_lin2

        print 'Best: SVC with RBF kernel'
        print best_rbf
        print combo_rbf

        print 'Best: SVC with polynomial (degree 4) kernel'
        print best_poly
        print combo_poly

        
        pl.show()



def amoc_aabw_timeseries_ext(ot):
        t_amoc = []; amoc = []; aabw = []; eq = []
        for i in range(len(ot.time)):
                t_amoc.append(ot.time[i].values/360.)
                amoc.append(ot.isel(time=i, depth=slice(0,32), lat=slice(24,40)).max().values/1000000.)
                eq.append(ot.isel(time=i, depth=slice(0,27), lat=22).max().values/1000000.)
                aabw.append(ot.isel(time=i, depth=slice(0,8), lat=slice(0,10)).min().values/1000000.)
        return t_amoc, amoc, eq, aabw

if __name__ == '__main__':
        MAIN()

