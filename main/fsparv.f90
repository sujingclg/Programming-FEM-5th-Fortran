SUBROUTINE fsparv(kv,km,g,kdiag)
!
! This subroutine assembles element matrices into a symmetric skyline
! global matrix.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(IN)::g(:),kdiag(:)
 REAL(iwp),INTENT(IN)::km(:,:)
 REAL(iwp),INTENT(OUT)::kv(:) 
 INTEGER::i,idof,k,j,iw,ival
 idof=UBOUND(g,1) ! 单元自由度
 DO i=1,idof
   k=g(i) ! 自由度编号
   IF(k/=0)THEN ! 说明是无约束
     DO j=1,idof
       IF(g(j)/=0)THEN
         iw=k-g(j)
         IF(iw>=0)THEN
           ival=kdiag(k)-iw
           kv(ival)=kv(ival)+km(i,j) ! 将总刚度矩阵压缩到一维数组上，每次迭代会在单个元素上面累加进行单刚组装
         END IF
       END IF
     END DO
   END IF
 END DO
RETURN
END SUBROUTINE fsparv