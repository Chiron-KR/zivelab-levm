# Battery full cell

전극과 전해액 사이의 인터페이스에서는 (1) charge transfer reaction, (2) 농도 차이에 의해 발생하는 diffusion, (3) double layer capacitor, (4) 전해액 저항까지 4가지 요소를 고려해야 합니다. 이를 나타내는 등가회로는 Randles circuit입니다. 

	Z = Rs - (Rct - W) | Cdl
	
위 표현식은 diffusion으로 Warburg impedance W를 사용하였다. 저 주파수에서 -45도의 위상을 가지고 Z', Z"이 같이 커진다.

## Cathode

Li-ion battery의 양극에 대한 half cell impedance를 측정하면 diffusion의 형태가 다릅니다. 저 주파수에서 마치 Capacitor인 것처럼 행동한다. 복잡한 설명은 버리고 그냥 이런 형태의 diffusion 모델을 T라 하자.  그러면 이 half cell의 임피던스는 아래처럼 쓸 수 있다.  

	Z = Rs - (Rct - T) | Cdl
	
![Cathode](https://raw.githubusercontent.com/zivelab/zivelab-levm/master/doc/images/Cathode.png)
Rs = 80mΩ, Rct = 1Ω, Ty = 1/sqrt(Rd/Cd) and Tb = sqrt(Rd*Cd), Cdl = 1mΩ where Rd = 7Ω/cm and Cd = 2F/cm 

## Anode

반면, Li metal 같은 음극에 대한 half cell에서는 고주파수에서 W(-45도 위상)처럼 행동하고, 저 주파수에서 R | C (반원)처럼 행동하는 임피던스를 측정하게 됩니다. 이 때 diffusion 모델을 O라 하자. 마찬가지로 임피던스는 아래처럼 쓸 수 있다. 

	Z = Rs - (Rct - O) | Cdl
	
![Anode](https://raw.githubusercontent.com/zivelab/zivelab-levm/master/doc/images/Anode.png)
Rs = 80mΩ, Rct = 1Ω, Oy = 1/sqrt(rd/Cd) and Ob = sqrt(Rd*Cd), Cdl = 10mΩ where Rd = 7Ω/cm and Cd = 1F/cm

## Diffusion Elements

아래 그림들은 W, T, O에 대한 전형적인 모습을 보여준다.
	
### W element

![W](https://raw.githubusercontent.com/zivelab/zivelab-levm/master/doc/images/W.png)
W = 1/sqrt(Rd/Cd) where Rd = 1Ω/cm, Cd = 1F/cm

### T element
![T](https://raw.githubusercontent.com/zivelab/zivelab-levm/master/doc/images/T.png)
Ty = 1/sqrt(Rd/Cd), Tb = sqrt(Rd*Cd) where Rd = 1Ω, Cd = 300F

###  O element
![O](https://raw.githubusercontent.com/zivelab/zivelab-levm/master/doc/images/O.png)
Oy = 1/sqrt(Rd/Cd), Ob = sqrt(Rd*Cd) where Rd = 1Ω, Cd = 10F

## Battery Full Cell

그럼 full cell에 대한 임피던스는 전해액 저항 + 양극 저항 + 음극 저항이 되므로 다음과 같이 쓸 수 있다. 

	Z = Ls - Rs - (Rct_c - T_c) | Cdl_c - (Rct_a - O_a) | Cdl_a
	
여기서 각 파라미터의 물리적 의미는 다음과 같다.

- `Ls` : geometric inductivity
- `Rs` : Electrolyte resistance
- `Rct_c` : Charge transfer resistance of cathode
- `Cdl_c` : Double layer capacitance of cathode
- `T_c` : finite length reflective diffusion inside cathode particles
- `Rct_a` : Charge transfer resistance of anode
- `Cdl_a` : Double layer capacitance of anode
- `O_a` : finite length transmissive diffusion through passivating layer on anode surface

전형적인 모습은 아래와 같다.

![full cell](https://raw.githubusercontent.com/zivelab/zivelab-levm/master/doc/images/FullCell.png)
Ls = 5uH, Rs = 0.04Ω, Rct_c = 0.4Ω, Ty_c = 25.8, Tb_c = 77.46, Cdl_c = 0.01F, Rct_a = 0.2Ω, Oy_a = 44.7, Ob_a = 22.36, Cdl_a = 1mF

## Simplified Battery Full Cell

T와 O는 Hyperbolic Cotangent, Tangent 함수가 필요한 계산이 비싼 함수이다. 그던데, 저주파수에서 각각 C, R | C로 근사될 수 있으므로 아래와 같이 계산하기 쉬운 수학식으로 쓸 수 있다.

	Z = Ls - Rs - (Rct_c - Cd_c) | Cdl_c - (Rct_a - Rd_a | Cd_a) | Cdl_a
	
- `C_c` : total intercalation capacitance of cathode
- `Rd_a` : low frequency limit of finite length transmissive diffusion through passivating layer on anode surface
- `Cd_a` : capacitance due to dielectric relaxation of passivating layer

![simplified full cell](https://raw.githubusercontent.com/zivelab/zivelab-levm/master/doc/images/SimpleFullCell.png)
Cd_c=2000F, Rd_a=0.5Ω, Cd_a=1000F

	







	

	
	