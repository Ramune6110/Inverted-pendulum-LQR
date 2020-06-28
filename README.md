# Inverted-pendulum-LQR
Inverted pendulum LQR
## Constitution
lqr_control.hにclassを定義.
関数の実装はlqr_control.cppに書いてあります.
またmain.cppがmainのファイルになります.
実行方法は以下のProcedureを参考に行ってください.
実行にはEigenのインストールが必要です.
```bash
sudo apt-get install libeigen3-dev
```
## Result
![LQR](https://github.com/Ramune6110/Inverted-pendulum-LQR/blob/master/LQR_result.png)  
## Environment
Ubuntu18.04
## Procedure
```bash
g++ main.cpp lqr_control.cpp -I /usr/include/eigen3
```
```bash
./a.out
```

