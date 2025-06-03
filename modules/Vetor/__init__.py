

import math
import matplotlib.pyplot as plt


class Vetor:
    def __init__(self, dados):
        if not isinstance(dados, list):
            raise TypeError("A entrada deve ser uma lista.")
        
        if not all(isinstance(x, (int, float)) for x in dados):
            raise ValueError("Todos os elementos da lista devem ser int ou float.")
        
        if not 1 <= len(dados) <= 3:
            raise ValueError("A lista deve conter entre 1 e 3 elementos.")
        
        self.dim = len(dados)
        self.x = dados[0]
        self.y = dados[1] if len(dados) > 1 else None
        self.z = dados[2] if len(dados) > 2 else None

    def __str__(self):
        componentes = [self.x]
        if self.y is not None:
            componentes.append(self.y)
        if self.z is not None:
            componentes.append(self.z)
        return f"Vetor{tuple(componentes)}"

    def __add__(self, outro):
            if not isinstance(outro, Vetor):
                return NotImplemented
            if self.dim != outro.dim:
                raise ValueError("Vetores com dimensões diferentes não podem ser somados.")
            if self.dim == 1:
                return Vetor([self.x + outro.x])
            elif self.dim == 2:
                return Vetor([self.x + outro.x, self.y + outro.y])
            else:
                return Vetor([self.x + outro.x, self.y + outro.y, self.z + outro.z])

    def __mul__(self, escalar):
        if not isinstance(escalar, (int, float)):
            return NotImplemented
        if self.dim == 1:
            return Vetor([self.x * escalar])
        elif self.dim == 2:
            return Vetor([self.x * escalar, self.y * escalar])
        else:
            return Vetor([self.x * escalar, self.y * escalar, self.z * escalar])

    def __rmul__(self, escalar):
        return self.__mul__(escalar)

    def __sub__(self, outro):
            if not isinstance(outro, Vetor):
                return NotImplemented
            if self.dim != outro.dim:
                raise ValueError("Vetores com dimensões diferentes não podem ser subtraídos.")
            if self.dim == 1:
                return Vetor([self.x - outro.x])
            elif self.dim == 2:
                return Vetor([self.x - outro.x, self.y - outro.y])
            else:
                return Vetor([self.x - outro.x, self.y - outro.y, self.z - outro.z])

    def __neg__(self):
        if self.dim == 1:
            return Vetor([-self.x])
        elif self.dim == 2:
            return Vetor([-self.x, -self.y])
        else:
            return Vetor([-self.x, -self.y, -self.z])


    def modulo(self):
        soma = 0
        soma += self.x * self.x if self.x is not None else 0
        soma += self.y * self.y if self.y is not None else 0
        soma += self.z * self.z if self.z is not None else 0
        return math.sqrt(soma)

    

    def produto_interno(self, outro):
        if not isinstance(outro, Vetor):
            return NotImplemented
        if self.dim != outro.dim:
            raise ValueError("Vetores com dimensões diferentes.")
        soma = 0
        soma += self.x * outro.x if self.x is not None else 0
        soma += self.y * outro.y if self.y is not None else 0
        soma += self.z * outro.z if self.z is not None else 0
        return soma

    def angulo(self, outro):
        if not isinstance(outro, Vetor):
            return NotImplemented
        if self.dim != outro.dim:
            raise ValueError("Vetores com dimensões diferentes.")
        mod_self = self.modulo()
        mod_outro = outro.modulo()
        if mod_self == 0 or mod_outro == 0:
            raise ValueError("Vetor Nulo")
        cos_theta = self.produto_interno(outro) / (mod_self * mod_outro)
        angulo_rad = math.acos(cos_theta)
        return math.degrees(angulo_rad)

    def area_ao_quadrado(self, outro):
        area = self.produto_interno(self) * outro.produto_interno(outro) - self.produto_interno(outro) ** 2
        return area

    def det_2d(self, outro):
        return self.x * outro.y - self.y * outro.x

    def projecao(self, outro):
        return (self.produto_interno(outro) / self.produto_interno(self)) * self

    def produto_vetorial(self, outro):
        if not isinstance(outro, Vetor):
            return NotImplemented
        if self.dim != 3 or outro.dim != 3:
            raise ValueError("O produto vetorial só é definido para vetores de dimensão 3.")
        
        cx = self.y * outro.z - self.z * outro.y
        cy = self.z * outro.x - self.x * outro.z
        cz = self.x * outro.y - self.y * outro.x

        return Vetor([cx, cy, cz])

    def det_3d(self, v2, v3):
        if not all(isinstance(v, Vetor) for v in (v2, v3)):
            raise TypeError("Todos os argumentos devem ser Vetor.")
        if not all(v.dim == 3 for v in (self, v2, v3)):
            raise ValueError("O determinante 3D só é definido para vetores de dimensão 3.")
        
        a1, a2, a3 = self.x, self.y, self.z
        b1, b2, b3 = v2.x, v2.y, v2.z
        c1, c2, c3 = v3.x, v3.y, v3.z

        return (
            a1 * (b2 * c3 - b3 * c2)
            - a2 * (b1 * c3 - b3 * c1)
            + a3 * (b1 * c2 - b2 * c1)
        )

    def plot(self, fig=None, cor='b', limite=1.5, mostrar_eixos=True):
        if self.dim == 1:
            raise NotImplementedError("Plot não implementado para vetores de dimensão 1.")
        
        if fig is None:
            fig = plt.figure()
            if self.dim == 2:
                ax = fig.add_subplot(1,1,1)
            else:
                ax = fig.add_subplot(1,1,1, projection='3d')
        else:
            ax = fig.gca()

        if self.dim == 2:
            ax.quiver(0, 0, self.x, self.y, angles='xy', scale_units='xy', scale=1, color=cor)

            if mostrar_eixos:
                max_range = max(abs(self.x), abs(self.y)) * limite + 1e-6
                ax.axhline(0, color='k', linewidth=0.5)
                ax.axvline(0, color='k', linewidth=0.5)
                ax.set_xlim(-max_range, max_range)
                ax.set_ylim(-max_range, max_range)
                ax.set_aspect('equal')

        else:
            ax.quiver(0, 0, 0, self.x, self.y, self.z, color=cor)

            if mostrar_eixos:
                max_range = max(abs(self.x), abs(self.y), abs(self.z)) * limite + 1e-6
                ax.plot([-max_range, max_range], [0,0], [0,0], color='k', linewidth=0.5)
                ax.plot([0,0], [-max_range, max_range], [0,0], color='k', linewidth=0.5)
                ax.plot([0,0], [0,0], [-max_range, max_range], color='k', linewidth=0.5)

                ax.set_xlim([-max_range, max_range])
                ax.set_ylim([-max_range, max_range])
                ax.set_zlim([-max_range, max_range])

        return fig

    def distancia(self, outro):
        from modules.Reta import Reta
        if isinstance(outro, Vetor):
            if self.dim != outro.dim:
                raise ValueError("Vetores com dimensões diferentes.")
            return (self - outro).modulo()
        
        elif isinstance(outro, Reta):
            if self.dim != outro.dim:
                raise ValueError("A dimensão do vetor e da reta deve coincidir.")

            P = self
            A = outro.ponto
            v = outro.diretor
            n = outro.normal


            if self.dim == 2:
                # Sempre existe vetor normal
                a, b = n.x, n.y
                c = n.produto_interno(A)
                numerador = abs(a * self.x + b * self.y - c)
                denominador = math.sqrt(a**2 + b**2)
                return numerador / denominador

            elif self.dim == 3:
                # Distância ponto-reta em 3D: ||(P - A) × V|| / ||V||
                vetor = P - A
                prod_vetorial = vetor.produto_vetorial(v)
                return prod_vetorial.modulo() / v.modulo()
        
        else:
            raise TypeError("distancia_ate espera um objeto Vetor ou Reta.")


    # def distancia(self, outro):
